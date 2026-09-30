import os
import json
import time
import requests
import pandas as pd
import gspread
from google.oauth2.service_account import Credentials

# 根據 ASM mBio (10.1128/mbio.02013-25) 整理之主要 RNA 病毒清單
ORGANISM_LIST = [
    # Coronaviridae (冠狀病毒科)
    'SARS-CoV-2', 
    'MERS-CoV', 
    'SARS-CoV-1', 
    'Human coronavirus 229E', 
    'Human coronavirus NL63', 
    'Human coronavirus OC43',
    'Human coronavirus HKU1',

    # Flaviviridae (黃病毒科)
    'Hepatitis C virus', 
    'Dengue virus', 
    'Yellow fever virus', 
    'Zika virus', 
    'West Nile virus', 
    'Japanese encephalitis virus',

    # Picornaviridae (微小RNA病毒科)
    'Hepatitis A virus', 
    'Coxsackievirus', 
    'Enterovirus A71',
    'Rhinovirus',

    # Caliciviridae & Togaviridae (杯狀與披膜病毒科)
    'Norovirus', 
    'Chikungunya virus', 
    'Rubella virus',

    # Orthomyxoviridae & Paramyxoviridae (正黏與副黏液病毒科)
    'Influenza A virus', 
    'Influenza B virus',
    'Respiratory syncytial virus', 
    'Mumps virus', 
    'Measles virus', 
    'Human metapneumovirus',

    # Rhabdoviridae, Filoviridae & Arenaviridae (線狀、絲狀與沙狀病毒科)
    'Rabies virus', 
    'Ebola virus', 
    'Marburg virus',
    'Lassa virus', 
    'Crimean-Congo hemorrhagic fever virus',

    # Reoviridae (重組RNA病毒科)
    'Rotavirus'
]

SCOPES = [
    'https://www.googleapis.com/auth/spreadsheets',
    'https://www.googleapis.com/auth/drive'
]

# 計數器：用於實現「每 3 次查詢暫停 3 秒」
request_counter = 0

def rate_limit_control():
    """控制每查詢 3 次自動暫停 3 秒"""
    global request_counter
    request_counter += 1
    if request_counter % 3 == 0:
        print("  ⏳ [Rate Limit] 已執行 3 次查詢，自動暫停 3 秒...")
        time.sleep(3)

def get_cids_by_virus(virus_name, max_results=20):
    """根據 RNA 病毒名稱搜尋 PubChem 中對應的 RdRP 抑制劑 CID"""
    rate_limit_control()
    query = f"{virus_name} RdRP inhibitor"
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/cids/JSON"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            return data.get("IdentifierList", {}).get("CID", [])[:max_results]
    except Exception as e:
        print(f"[{virus_name}] 搜尋 CID 失敗: {e}")
    return []

def get_compound_structures(cid):
    """獲取化合物結構資訊：Canonical SMILES 與 InChIKey"""
    rate_limit_control()
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES,InChIKey/JSON"
    try:
        res = requests.get(url, timeout=10)
        if res.status_code == 200:
            props = res.json().get("PropertyTable", {}).get("Properties", [])[0]
            return props.get("CanonicalSMILES", ""), props.get("InChIKey", "")
    except Exception:
        pass
    return "", ""

def get_bioassay_data(cid):
    """取得 CID 的生物活性數據，僅留下明確有抑制劑數值 (IC50, EC50, Ki, Kd) 的紀錄"""
    rate_limit_control()
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    records = []
    try:
        res = requests.get(url, timeout=10)
        if res.status_code == 200:
            columns = res.json().get("Table", {}).get("ColumnName", [])
            rows = res.json().get("Table", {}).get("Row", [])
            
            def idx(col_name):
                return columns.index(col_name) if col_name in columns else -1

            i_target = idx("TargetName")
            i_type = idx("ActivityType")
            i_val = idx("ActivityValue")
            i_unit = idx("ActivityUnit")
            i_relation = idx("ActivityRelation")
            i_gene = idx("TargetGeneID")

            for row in rows:
                cell = row.get("Cell", [])
                activity_type = cell[i_type] if i_type != -1 and i_type < len(cell) else ""
                activity_val = cell[i_val] if i_val != -1 and i_val < len(cell) else ""
                
                # 只有活性數值不為空，且屬於 IC50/EC50/Ki/Kd 才視為有效抑制劑紀錄
                if activity_type.upper() in ["IC50", "EC50", "KI", "KD"] and str(activity_val).strip() != "":
                    records.append({
                        "Target Name": cell[i_target] if i_target != -1 and i_target < len(cell) else "",
                        "Standard Type": activity_type,
                        "Standard Value": activity_val,
                        "Standard Units": cell[i_unit] if i_unit != -1 and i_unit < len(cell) else "",
                        "Standard Relation": cell[i_relation] if i_relation != -1 and i_relation < len(cell) else "",
                        "NCBI Gene ID": cell[i_gene] if i_gene != -1 and i_gene < len(cell) else ""
                    })
    except Exception as e:
        print(f"獲取 CID {cid} Bioassay 失敗: {e}")
    return records

def get_uniprot_and_protein_id(gene_id, virus_name):
    """交叉比對 NCBI Gene ID 與 UniProt ID"""
    if gene_id and gene_id != "N/A" and str(gene_id).strip() != "":
        rate_limit_control()
        url = f"https://rest.uniprot.org/uniprotkb/search?query=xref:geneid-{gene_id}&fields=accession,id,protein_name"
        try:
            res = requests.get(url, timeout=10)
            if res.status_code == 200:
                results = res.json().get("results", [])
                if results:
                    return results[0].get("primaryAccession", ""), f"NCBI_Gene:{gene_id}"
        except Exception:
            pass
        return "", f"NCBI_Gene:{gene_id}"

    rate_limit_control()
    clean_virus = virus_name.replace("–", "-")
    url = f"https://rest.uniprot.org/uniprotkb/search?query=organism_name:\"{clean_virus}\"+AND+rdrp&fields=accession"
    try:
        res = requests.get(url, timeout=10)
        if res.status_code == 200:
            results = res.json().get("results", [])
            if results:
                return results[0].get("primaryAccession", ""), "N/A"
    except Exception:
        pass

    return "N/A", "N/A"

def update_google_sheet(df):
    """按病毒名稱分組，僅將「確定有抑制劑數據」的病毒寫入 Google Sheet 獨立工作表」"""
    if df.empty:
        print("⚠️ 未找到任何帶有有效活性數據的抑制劑，取消更新 Google Sheet。")
        return

    creds_json_str = os.environ.get("GCP_SA_KEY")
    spreadsheet_id = os.environ.get("GOOGLE_SHEET_ID")

    if not creds_json_str or not spreadsheet_id:
        raise ValueError("缺失 GCP_SA_KEY 或 GOOGLE_SHEET_ID 環境變數！")

    creds_dict = json.loads(creds_json_str)
    creds = Credentials.from_service_account_info(creds_dict, scopes=SCOPES)
    gc = gspread.authorize(creds)

    sh = gc.open_by_key(spreadsheet_id)
    existing_worksheets = {ws.title: ws for ws in sh.worksheets()}

    for organism, group_df in df.groupby("Organism"):
        sheet_title = str(organism).replace("–", "-")[:31]
        
        if sheet_title in existing_worksheets:
            worksheet = existing_worksheets[sheet_title]
        else:
            worksheet = sh.add_worksheet(title=sheet_title, rows=100, cols=20)
            existing_worksheets[sheet_title] = worksheet
        
        worksheet.clear()

        df_clean = group_df.fillna("")
        data_to_write = [df_clean.columns.values.tolist()] + df_clean.values.tolist()
        
        worksheet.update(data_to_write)
        print(f"✅ 已成功將 {len(group_df)} 筆抑制劑數據存入工作表 [{sheet_title}]")

    print(f"\n🎉 含有抑制劑數據的病毒分頁更新完畢！(試算表 ID: {spreadsheet_id})")

def main():
    all_rows = []
    
    for virus in ORGANISM_LIST:
        print(f"\n正在查詢 RNA 病毒: {virus} ...")
        cids = get_cids_by_virus(virus)
        
        virus_inhibitors_found = 0
        
        for cid in cids:
            assays = get_bioassay_data(cid)
            
            # 關鍵修改：只有當 assays 非空（代表該 CID 確實有發布抑制劑活性數據）時才進行處理解析
            if assays:
                smiles, inchikey = get_compound_structures(cid)
                for assay in assays[:2]: # 每個化合物最多擷取前 2 筆關聯性高的數據
                    gene_id = assay["NCBI Gene ID"]
                    uniprot_id, ncbi_protein_id = get_uniprot_and_protein_id(gene_id, virus)
                    
                    all_rows.append({
                        "Organism": virus,
                        "Target Name": assay["Target Name"] or "RdRP Polymerase",
                        "UniProt ID": uniprot_id,
                        "NCBI Gene ID": gene_id or "N/A",
                        "NCBI Protein ID": ncbi_protein_id,
                        "PubChem CID": str(cid),
                        "Standard Type": assay["Standard Type"],
                        "Standard Value": str(assay["Standard Value"]),
                        "Standard Units": assay["Standard Units"],
                        "Standard Relation": assay["Standard Relation"],
                        "Canonical SMILES": smiles,
                        "InChIKey": inchikey
                    })
                    virus_inhibitors_found += 1

        if virus_inhibitors_found > 0:
            print(f"  👉 [{virus}] 找到 {virus_inhibitors_found} 筆有效抑制劑資料，準備寫入！")
        else:
            print(f"  ⚪ [{virus}] 未發現帶有明確活性數據的抑制劑，跳過此病毒。")

    df = pd.DataFrame(all_rows)
    update_google_sheet(df)

if __name__ == "__main__":
    main()
