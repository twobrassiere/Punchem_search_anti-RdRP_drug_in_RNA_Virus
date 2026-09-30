import os
import json
import time
import random
import requests
import pandas as pd
import gspread
from google.oauth2.service_account import Credentials

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

request_counter = 0

def rate_limit_control():
    """控制每 3 次查詢自動暫停 3 秒"""
    global request_counter
    request_counter += 1
    if request_counter % 3 == 0:
        print("  ⏳ [Rate Limit] 已執行 3 次 API 查詢，自動暫停 3 秒...")
        time.sleep(3)

def fetch_ncbi_virus_full_name(virus_name):
    """【備援機制 1】向 NCBI Taxonomy / Entrez API 查詢病毒官方完整名稱與 TaxID"""
    rate_limit_control()
    print(f"  🔍 [NCBI Virus] 正在至 NCBI Taxonomy 檢索病毒完整名稱: {virus_name} ...")
    
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={virus_name}&retmode=json"
    try:
        res = requests.get(search_url, timeout=10)
        if res.status_code == 200:
            id_list = res.json().get("esearchresult", {}).get("idlist", [])
            if id_list:
                tax_id = id_list[0]
                rate_limit_control()
                summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id={tax_id}&retmode=json"
                sum_res = requests.get(summary_url, timeout=10)
                if sum_res.status_code == 200:
                    result = sum_res.json().get("result", {}).get(tax_id, {})
                    official_name = result.get("scientificname", virus_name)
                    print(f"  ✨ [NCBI Virus] 取得 NCBI 官方全名: {official_name} (TaxID: {tax_id})")
                    return official_name, tax_id
    except Exception as e:
        print(f"  ⚠️ [NCBI Virus] 查詢失敗: {e}")
    
    return virus_name, "N/A"

def get_cids_by_virus(query_term, max_results=20):
    """根據關鍵字搜尋 PubChem CID"""
    rate_limit_control()
    query = f"{query_term} RdRP inhibitor"
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/cids/JSON"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            return data.get("IdentifierList", {}).get("CID", [])[:max_results]
    except Exception:
        pass
    return []

def get_compound_structures(cid):
    """獲取化合物結構資訊：SMILES 與 InChIKey"""
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
    """取得 CID 在 PubChem 的生物活性數據 (IC50, EC50, Ki, Kd)"""
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
                
                if activity_type.upper() in ["IC50", "EC50", "KI", "KD"] and str(activity_val).strip() != "":
                    records.append({
                        "Target Name": cell[i_target] if i_target != -1 and i_target < len(cell) else "",
                        "Standard Type": activity_type,
                        "Standard Value": activity_val,
                        "Standard Units": cell[i_unit] if i_unit != -1 and i_unit < len(cell) else "",
                        "Standard Relation": cell[i_relation] if i_relation != -1 and i_relation < len(cell) else "",
                        "NCBI Gene ID": cell[i_gene] if i_gene != -1 and i_gene < len(cell) else ""
                    })
    except Exception:
        pass
    return records

def get_chembl_activity_data(virus_name, max_results=5):
    """透過 ChEMBL API 檢索該病毒 RdRP 的藥物活性數據"""
    rate_limit_control()
    print(f"  🧪 [ChEMBL] 正在檢索 ChEMBL 資料庫: {virus_name} ...")
    records = []
    
    clean_virus = virus_name.replace("–", "-")
    url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_organism__icontains={clean_virus}&standard_type__in=IC50,EC50,Ki,Kd&limit={max_results}"
    
    try:
        res = requests.get(url, timeout=12)
        if res.status_code == 200:
            activities = res.json().get("activities", [])
            for act in activities:
                molecule_chembl_id = act.get("molecule_chembl_id", "N/A")
                standard_type = act.get("standard_type", "")
                standard_value = act.get("standard_value", "")
                standard_units = act.get("standard_units", "")
                standard_relation = act.get("standard_relation", "=")
                target_pref_name = act.get("target_pref_name", "RdRP Polymerase")
                canonical_smiles = act.get("canonical_smiles", "N/A")
                
                if standard_value and standard_type:
                    records.append({
                        "Organism": virus_name,
                        "Target Name": f"{target_pref_name} (ChEMBL)",
                        "UniProt ID": act.get("target_chembl_id", "N/A"),
                        "NCBI Gene ID": "N/A",
                        "NCBI Protein ID": "N/A",
                        "PubChem CID": f"ChEMBL:{molecule_chembl_id}",
                        "Standard Type": standard_type,
                        "Standard Value": str(standard_value),
                        "Standard Units": standard_units or "nM",
                        "Standard Relation": standard_relation,
                        "Canonical SMILES": canonical_smiles or "N/A",
                        "InChIKey": "N/A"
                    })
            if records:
                print(f"  ✨ [ChEMBL] 成功取得 {len(records)} 筆 ChEMBL 活性數據！")
    except Exception as e:
        print(f"  ⚠️ [ChEMBL] 檢索失敗: {e}")
        
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
    """將單一病毒的最新數據更新至 Google Sheet 中的專屬工作表 (Worksheet)"""
    if df.empty:
        print("⚠️ 未抓取到任何資料，取消更新 Google Sheet。")
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
        
        # 清空該病毒對應的工作表並寫入最新搜尋結果
        worksheet.clear()

        df_clean = group_df.fillna("")
        data_to_write = [df_clean.columns.values.tolist()] + df_clean.values.tolist()
        
        worksheet.update(data_to_write)
        print(f"✅ 已成功更新工作表 [{sheet_title}] (共 {len(group_df)} 筆資料)")

    print(f"\n🎉 Google Sheet 單一病毒工作表同步完畢！(試算表 ID: {spreadsheet_id})")

def generate_email_summary(df, virus_name):
    """分析本次抽查的病毒結果，並在 Console Log 輸出總結"""
    total_records = len(df)
    active_drugs = df[df["PubChem CID"] != "No Active Inhibitor Found"]
    
    report_lines = [
        "==========================================",
        "  PubChem / ChEMBL 隨機抽查病毒每日報告",
        "==========================================",
        f"🎲 本次抽查病毒: {virus_name}",
        f"📅 執行時間: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}",
        f"📊 總檢索紀錄數: {total_records} 筆",
        f"💊 有效化合物數量: {active_drugs['PubChem CID'].nunique()} 個",
        "------------------------------------------"
    ]
    
    if not active_drugs.empty:
        report_lines.append("【活性類型 (Standard Type) 分布】:")
        type_counts = active_drugs["Standard Type"].value_counts()
        for stype, count in type_counts.items():
            report_lines.append(f"  • {stype}: {count} 筆試驗數據")

    report_lines.extend([
        "------------------------------------------",
        "🔗 最新數據已同步至 Google Sheet 對應的工作表。",
        "=========================================="
    ])
    
    report_text = "\n".join(report_lines)
    print("\n" + report_text)
    return report_text

def main():
    # 🎲 關鍵改動：每次執行時，從 ORGANISM_LIST 隨機挑選 1 種病毒
    selected_virus = random.choice(ORGANISM_LIST)
    
    print(f"==========================================")
    print(f"🎲 隨機抽中本次檢索病毒: {selected_virus}")
    print(f"==========================================")
    
    all_rows = []
    
    # 1. 第一階段：以原本名稱查詢 PubChem
    cids = get_cids_by_virus(selected_virus)
    assays_found = []

    for cid in cids:
        assays = get_bioassay_data(cid)
        if assays:
            smiles, inchikey = get_compound_structures(cid)
            for assay in assays[:2]:
                gene_id = assay["NCBI Gene ID"]
                uniprot_id, ncbi_protein_id = get_uniprot_and_protein_id(gene_id, selected_virus)
                assays_found.append({
                    "Organism": selected_virus,
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

    # 2. 第二階段：併行/備援查詢 ChEMBL 資料庫
    chembl_records = get_chembl_activity_data(selected_virus)
    if chembl_records:
        assays_found.extend(chembl_records)

    # 3. 第三階段：若 PubChem 與 ChEMBL 皆無活性藥物，觸發 NCBI Virus (Entrez API) 檢索
    if not assays_found:
        print(f"  ⚪ [PubChem/ChEMBL] 原名稱無活性藥物，啟動 NCBI Virus 備援機制...")
        ncbi_official_name, tax_id = fetch_ncbi_virus_full_name(selected_virus)
        
        if ncbi_official_name != selected_virus:
            retry_cids = get_cids_by_virus(ncbi_official_name)
            for cid in retry_cids:
                retry_assays = get_bioassay_data(cid)
                if retry_assays:
                    smiles, inchikey = get_compound_structures(cid)
                    for assay in retry_assays[:2]:
                        gene_id = assay["NCBI Gene ID"]
                        uniprot_id, ncbi_protein_id = get_uniprot_and_protein_id(gene_id, ncbi_official_name)
                        assays_found.append({
                            "Organism": selected_virus,
                            "Target Name": f"{assay['Target Name']} (NCBI: {ncbi_official_name})",
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

        if not assays_found:
            uniprot_id, _ = get_uniprot_and_protein_id("", ncbi_official_name)
            assays_found.append({
                "Organism": selected_virus,
                "Target Name": f"RNA-dependent RNA polymerase ({ncbi_official_name})",
                "UniProt ID": uniprot_id,
                "NCBI Gene ID": f"NCBI_TaxID:{tax_id}",
                "NCBI Protein ID": "N/A",
                "PubChem CID": "No Active Inhibitor Found",
                "Standard Type": "N/A",
                "Standard Value": "N/A",
                "Standard Units": "N/A",
                "Standard Relation": "N/A",
                "Canonical SMILES": "N/A",
                "InChIKey": "N/A"
            })

    all_rows.extend(assays_found)

    df = pd.DataFrame(all_rows)
    update_google_sheet(df)
    generate_email_summary(df, selected_virus)

if __name__ == "__main__":
    main()
