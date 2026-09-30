import os
import json
import time
import re
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
    """配合 API 規範微幅延遲，避免 API 請求過快失敗"""
    global request_counter
    request_counter += 1
    if request_counter % 3 == 0:
        time.sleep(0.5)

def extract_mutations_from_text(text):
    """從文獻全文/摘要中提取 RdRP 胺基酸突變點 (如 S759A, V557L, F480L)"""
    if not text:
        return "None Found"
    pattern = r'\b[A-Z]\d{2,4}[A-Z]\b'
    matches = re.findall(pattern, text)
    filtered = [m for m in set(matches) if not m.startswith(('IC', 'EC', 'KI', 'KD', 'DNA', 'RNA', 'PCR'))]
    return ", ".join(filtered) if filtered else "None Detected"

def extract_bioactivity_from_text(text):
    """從文獻摘要中抓取 IC50 / EC50 生物活性數值與單位"""
    if not text:
        return "N/A", "N/A", "N/A", "="
    
    pattern = r'\b(IC50|EC50|Ki|Kd)\b\s*(=|:|of|is|around)?\s*([0-9\.]+\s*(?:uM|nM|mM|µM))'
    match = re.search(pattern, text, re.IGNORECASE)
    if match:
        act_type = match.group(1).upper()
        relation = match.group(2) if match.group(2) in ['=', ':'] else "="
        val_unit = match.group(3).strip().split()
        val = val_unit[0] if len(val_unit) > 0 else "N/A"
        unit = val_unit[1] if len(val_unit) > 1 else "uM"
        return act_type, val, unit, relation
    
    return "N/A", "N/A", "N/A", "="

def search_smiles_by_compound_name(compound_name):
    """利用化合物名稱查詢 PubChem 取得 Canonical SMILES 與 InChIKey"""
    if not compound_name or len(compound_name) < 3:
        return "N/A", "N/A"
    rate_limit_control()
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/CanonicalSMILES,InChIKey/JSON"
    try:
        res = requests.get(url, timeout=8)
        if res.status_code == 200:
            props = res.json().get("PropertyTable", {}).get("Properties", [])[0]
            return props.get("CanonicalSMILES", "N/A"), props.get("InChIKey", "N/A")
    except Exception:
        pass
    return "N/A", "N/A"

def fetch_pubmed_until_found(virus_name):
    """【無上限全量搜尋】持續翻頁查詢 PubMed API，直到掃描完所有文獻，搜集所有具備活性資訊或結構的資料為止"""
    print(f"  📚 [PubMed] 啟動全量檢索模式：掃描 {virus_name} RdRP inhibitor 所有文獻...")
    records = []
    
    term = f"({virus_name}[Title/Abstract]) AND (RdRP inhibitor)"
    retstart = 0
    batch_size = 20  # 每頁拉取 20 篇進行 XML 解析
    total_found_in_pubmed = None

    while True:
        rate_limit_control()
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={term}&retmode=json&retmax={batch_size}&retstart={retstart}"
        
        try:
            res = requests.get(search_url, timeout=10)
            if res.status_code != 200:
                break
                
            esearch_data = res.json().get("esearchresult", {})
            if total_found_in_pubmed is None:
                total_found_in_pubmed = int(esearch_data.get("count", 0))
                print(f"  ℹ️ [PubMed] 數據庫中符合條件之文獻總數: {total_found_in_pubmed} 篇")

            id_list = esearch_data.get("idlist", [])
            if not id_list or retstart >= total_found_in_pubmed:
                print("  🏁 [PubMed] 已完成所有頁面掃描。")
                break
            
            pmids = ",".join(id_list)
            rate_limit_control()
            fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmids}&retmode=xml"
            fetch_res = requests.get(fetch_url, timeout=15)
            
            if fetch_res.status_code == 200:
                xml_content = fetch_res.text
                articles = xml_content.split('<PubmedArticle>')
                
                for article_xml in articles[1:]:
                    pmid_match = re.search(r'<PMID.*?>(.*?)</PMID>', article_xml)
                    title_match = re.search(r'<ArticleTitle>(.*?)</ArticleTitle>', article_xml, re.DOTALL)
                    abstract_match = re.search(r'<AbstractText.*?>(.*?)</AbstractText>', article_xml, re.DOTALL)
                    
                    if pmid_match and title_match:
                        pmid = pmid_match.group(1).strip()
                        title = title_match.group(1).strip()
                        abstract = abstract_match.group(1).strip() if abstract_match else ""
                        full_text = f"{title} {abstract}"
                        
                        # 1. 抓取生物活性資訊 (IC50, EC50, Ki, Kd)
                        act_type, val, unit, relation = extract_bioactivity_from_text(full_text)
                        
                        # 2. 抓取突變點 (RdRP Mutation)
                        mutations = extract_mutations_from_text(full_text)
                        
                        # 3. 提取標題/摘要中的化合物詞彙並獲取 SMILES 結構
                        smiles, inchikey = "N/A", "N/A"
                        words = re.findall(r'\b[A-Za-z0-9\-]{4,25}\b', full_text)
                        for word in words:
                            if word.lower() not in ['virus', 'rdrp', 'inhibitor', 'activity', 'antiviral', 'sars', 'cov', 'protein', 'cell', 'study', 'assay']:
                                smiles, inchikey = search_smiles_by_compound_name(word)
                                if smiles != "N/A":
                                    break

                        # 🎯 條件判斷：僅留下「具備有效生物活性數值」或「成功對照出 SMILES 結構」的文獻紀錄
                        if val != "N/A" or smiles != "N/A":
                            records.append({
                                "Organism": virus_name,
                                "Target Name": f"Literature: {title[:70]}...",
                                "RdRP Mutation": mutations,  # 📍 Target Name 下一欄
                                "UniProt ID": "N/A",
                                "NCBI Gene ID": f"PubMed_PMID:{pmid}",
                                "NCBI Protein ID": "N/A",
                                "PubChem CID": f"PubMed:{pmid}",
                                "Standard Type": act_type if act_type != "N/A" else "Literature Record",
                                "Standard Value": str(val),
                                "Standard Units": unit if unit != "N/A" else "PubMed Record",
                                "Standard Relation": relation,
                                "Canonical SMILES": smiles,
                                "InChIKey": inchikey
                            })

            retstart += batch_size  # 自動推進至下一頁，直到掃描完 count 總數

        except Exception as e:
            print(f"  ⚠️ [PubMed] 搜尋過程發生異常: {e}")
            break

    print(f"  ✨ [PubMed] 掃描完畢！共採集到 {len(records)} 筆具備活性資訊或結構的有效文獻資料！")
    return records

def get_cids_by_virus(query_term, max_results=20):
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
                        "RdRP Mutation": "N/A",  # 📍 Target Name 下一欄
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
                        "RdRP Mutation": "N/A",  # 📍 Target Name 下一欄
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
            worksheet = sh.add_worksheet(title=sheet_title, rows=100, cols=25)
            existing_worksheets[sheet_title] = worksheet
        
        worksheet.clear()

        df_clean = group_df.fillna("")
        data_to_write = [df_clean.columns.values.tolist()] + df_clean.values.tolist()
        
        worksheet.update(data_to_write)
        print(f"✅ 已成功更新工作表 [{sheet_title}] (共 {len(group_df)} 筆資料)")

    print(f"\n🎉 單選病毒工作表同步完畢！(試算表 ID: {spreadsheet_id})")

def generate_email_summary(df, selected_virus):
    total_records = len(df)
    active_drugs = df[df["PubChem CID"] != "No Active Inhibitor Found"]
    
    report_lines = [
        "==========================================",
        "  PubChem / ChEMBL / PubMed 全量搜尋每日分析報告",
        "==========================================",
        f"📅 執行時間: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}",
        f"🎯 今日抽樣病毒: {selected_virus}",
        f"📊 總檢索紀錄數: {total_records} 筆",
        "------------------------------------------\n",
        "【數據來源與活性類型 (Standard Type) 分布】:"
    ]
    
    type_counts = active_drugs["Standard Type"].value_counts()
    for stype, count in type_counts.items():
        report_lines.append(f"  • {stype}: {count} 筆紀錄")

    report_lines.extend([
        "\n------------------------------------------",
        "🔗 最新數據已完整同步至 Google Sheet 試算表。",
        "=========================================="
    ])
    
    report_text = "\n".join(report_lines)
    print("\n" + report_text)
    return report_text

def main():
    all_rows = []
    
    selected_virus = random.choice(ORGANISM_LIST)
    
    print(f"\n==========================================")
    print(f"🎲 今日隨機抽樣搜尋 RNA 病毒: {selected_virus} ...")
    print(f"==========================================")
    
    # 1. 第一階段：PubChem BioAssay
    cids = get_cids_by_virus(selected_virus)
    assays_found = []

    for cid in cids:
        assays = get_bioassay_data(cid)
        if assays:
            smiles, inchikey = get_compound_structures(cid)
            for assay in assays[:3]:
                gene_id = assay["NCBI Gene ID"]
                uniprot_id, ncbi_protein_id = get_uniprot_and_protein_id(gene_id, selected_virus)
                assays_found.append({
                    "Organism": selected_virus,
                    "Target Name": assay["Target Name"] or "RdRP Polymerase",
                    "RdRP Mutation": "N/A",  # 📍 放在 Target Name 下一欄
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

    # 2. 第二階段：ChEMBL
    chembl_records = get_chembl_activity_data(selected_virus)
    if chembl_records:
        assays_found.extend(chembl_records)

    # 3. 第三階段：若前兩者無資料，觸發 PubMed 全量翻頁檢索（無上限掃描完所有文獻，僅收錄具備活性或結構的資料）
    if not assays_found:
        print(f"  ⚪ [PubChem/ChEMBL] 無數據，開啟 PubMed 全量無上限連鎖檢索...")
        pubmed_records = fetch_pubmed_until_found(selected_virus)
        if pubmed_records:
            assays_found.extend(pubmed_records)

    all_rows.extend(assays_found)

    df = pd.DataFrame(all_rows)
    update_google_sheet(df)
    generate_email_summary(df, selected_virus)

if __name__ == "__main__":
    main()
