import os
import time
import logging
import random
import requests
import urllib.parse
import re
import pandas as pd
from datetime import datetime
from typing import List, Dict, Set, Optional, Tuple

import json
try:
    import gspread
    from google.oauth2.service_account import Credentials
    GSPREAD_AVAILABLE = True
except ImportError:
    GSPREAD_AVAILABLE = False

# 設定 Logging 紀錄
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)

# 參考 AntiviralDB 收錄之完整 RNA 病毒種類清單
ANTIVIRALDB_RNA_VIRUS_LIST = [
    "SARS-CoV-2",
    "SARS-CoV",
    "MERS-CoV",
    "Human coronavirus 229E",
    "Human coronavirus NL63",
    "Hepatitis C virus",
    "Dengue virus",
    "Zika virus",
    "West Nile virus",
    "Yellow fever virus",
    "Japanese encephalitis virus",
    "Enterovirus 71",
    "Poliovirus",
    "Rhinovirus",
    "Hepatitis A virus",
    "Coxsackievirus",
    "Norovirus",
    "Chikungunya virus",
    "Rubella virus",
    "Influenza A virus",
    "Influenza B virus",
    "Respiratory syncytial virus",
    "Nipah virus",
    "Measles virus",
    "Mumps virus",
    "Rabies virus",
    "Ebola virus",
    "Marburg virus",
    "Lassa virus",
    "Crimean-Congo hemorrhagic fever virus",
    "Rotavirus"
]

NCBI_ESEARCH_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


# =============================================================================
# PubChem 藥物註解器 Class (含 Parent/Component CID 追蹤與 SMILES 備援)
# =============================================================================
class PubChemDrugAnnotator:
    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def __init__(self, delay_between_batches: float = 0.25, max_retries: int = 3, timeout: int = 15):
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "PubChemDrugAnnotator/3.0 (Python requests client)",
            "Accept": "application/json"
        })
        self.delay_between_batches = delay_between_batches
        self.max_retries = max_retries
        self.timeout = timeout

    def _safe_get_json(self, url: str, params: Optional[dict] = None) -> Optional[dict]:
        for attempt in range(self.max_retries):
            try:
                response = self.session.get(url, params=params, timeout=self.timeout)
                if response.status_code == 404:
                    return None
                response.raise_for_status()
                return response.json()
            except Exception:
                if attempt < self.max_retries - 1:
                    time.sleep(2 ** (attempt + 1))
                else:
                    return None
        return None

    def get_smiles_from_cid(self, cid: int) -> Optional[str]:
        """從 CID 取得 SMILES，包含 Parent CID 與 CACTUS 備援"""
        if not cid:
            return None

        # 1. 直接查詢 SMILES
        url = f"{self.BASE_URL}/compound/cid/{cid}/property/CanonicalSMILES,IsomericSMILES,SMILES/JSON"
        data = self._safe_get_json(url)
        if data and 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
            props = data['PropertyTable']['Properties'][0]
            smiles = props.get('CanonicalSMILES') or props.get('IsomericSMILES') or props.get('SMILES')
            if smiles:
                return smiles

        # 2. 查詢 Parent CID
        parent_url = f"{self.BASE_URL}/compound/cid/{cid}/cids/JSON?cids_type=parent"
        parent_data = self._safe_get_json(parent_url)
        if parent_data and 'IdentifierList' in parent_data and 'CID' in parent_data['IdentifierList']:
            parent_cid = parent_data['IdentifierList']['CID'][0]
            if str(parent_cid) != str(cid):
                p_smiles_url = f"{self.BASE_URL}/compound/cid/{parent_cid}/property/CanonicalSMILES,IsomericSMILES/JSON"
                p_data = self._safe_get_json(p_smiles_url)
                if p_data and 'PropertyTable' in p_data and 'Properties' in p_data['PropertyTable']:
                    p_props = p_data['PropertyTable']['Properties'][0]
                    p_smiles = p_props.get('CanonicalSMILES') or p_props.get('IsomericSMILES')
                    if p_smiles:
                        return p_smiles

        # 3. NCI CIR 備援
        try:
            cactus_url = f"https://cactus.nci.nih.gov/chemical/structure/pubchem:{cid}/smiles"
            resp = self.session.get(cactus_url, timeout=5)
            if resp.status_code == 200 and resp.text:
                return resp.text.strip()
        except Exception:
            pass

        return None

    def fetch_compound_properties_batch(self, cids: List[int]) -> List[Dict]:
        if not cids:
            return []

        properties_to_fetch = "Title,IUPACName,MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChIKey"
        chunk_size = 100
        all_properties = []

        for i in range(0, len(cids), chunk_size):
            sub_cids = ",".join(map(str, cids[i:i + chunk_size]))
            url = f"{self.BASE_URL}/compound/cid/{sub_cids}/property/{properties_to_fetch}/JSON"

            data = self._safe_get_json(url)
            if data and "PropertyTable" in data and "Properties" in data["PropertyTable"]:
                props_list = data["PropertyTable"]["Properties"]
                for prop in props_list:
                    cid = prop.get("CID")
                    smiles = prop.get("CanonicalSMILES") or prop.get("IsomericSMILES")

                    if not smiles and cid:
                        smiles = self.get_smiles_from_cid(cid)

                    prop["Canonical SMILES by pubchem"] = smiles if smiles else ""
                    all_properties.append(prop)

            time.sleep(self.delay_between_batches)

        return all_properties


annotator = PubChemDrugAnnotator()


# =============================================================================
# NCBI Gene API 模組 (查詢病毒 RdRP Domain 基因)
# =============================================================================
def fetch_ncbi_gene_rdrp_info(virus_name: str) -> Dict[str, str]:
    """使用 NCBI Entrez API 檢索特定病毒中與 RdRP 相關的 Gene ID 與 Symbol"""
    term = f"{virus_name}[Organism] AND (RdRP OR RNA-dependent RNA polymerase)"
    search_params = {
        "db": "gene",
        "term": term,
        "retmode": "json",
        "retmax": 1
    }
    try:
        res = requests.get(NCBI_ESEARCH_BASE, params=search_params, timeout=10)
        if res.status_code == 200:
            id_list = res.json().get("esearchresult", {}).get("idlist", [])
            if id_list:
                gene_id = id_list[0]
                # 取得詳細 Gene 摘要
                sum_params = {"db": "gene", "id": gene_id, "retmode": "json"}
                sum_res = requests.get(NCBI_ESUMMARY_BASE, params=sum_params, timeout=10)
                if sum_res.status_code == 200:
                    summary_data = sum_res.json().get("result", {}).get(str(gene_id), {})
                    return {
                        "NCBI Gene ID": gene_id,
                        "NCBI Gene Symbol": summary_data.get("name", "N/A"),
                        "NCBI Gene Description": summary_data.get("description", "N/A")
                    }
    except Exception as e:
        logging.error(f"查詢 NCBI Gene 失敗 [{virus_name}]: {e}")
        
    return {
        "NCBI Gene ID": "N/A",
        "NCBI Gene Symbol": "N/A",
        "NCBI Gene Description": "N/A"
    }


def search_cids_via_esearch(query: str, retmax: int = 200) -> List[int]:
    params = {
        "db": "pccompound",
        "term": query,
        "retmode": "json",
        "retmax": retmax
    }
    try:
        response = requests.get(NCBI_ESEARCH_BASE, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            id_list = data.get("esearchresult", {}).get("idlist", [])
            return [int(cid) for cid in id_list]
    except Exception as e:
        logging.error(f"ESearch 失敗 [{query}]: {e}")
    return []


def search_cids_via_pug_name(term: str) -> List[int]:
    url = f"{PubChemDrugAnnotator.BASE_URL}/compound/name/{requests.utils.quote(term)}/cids/JSON"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            return data.get("IdentifierList", {}).get("CID", [])
    except Exception as e:
        logging.debug(f"PUG Name 搜尋無結果 [{term}]: {e}")
    return []


def update_google_sheet_by_virus(virus_df_dict: Dict[str, pd.DataFrame], sheet_name: str = "RNA_Virus_RdRP_Inhibitors") -> None:
    if not GSPREAD_AVAILABLE:
        logging.warning("未安裝 gspread，跳過 Google Sheet 上傳。")
        return

    creds_env = os.environ.get("GOOGLE_SERVICE_ACCOUNT_JSON")
    creds_filename = os.environ.get("GOOGLE_APPLICATION_CREDENTIALS", "credentials.json")
    scopes = ["https://www.googleapis.com/auth/spreadsheets", "https://www.googleapis.com/auth/drive"]

    try:
        if creds_env:
            creds = Credentials.from_service_account_info(json.loads(creds_env), scopes=scopes)
        elif os.path.exists(creds_filename):
            creds = Credentials.from_service_account_file(creds_filename, scopes=scopes)
        else:
            logging.warning("找不到憑證，跳過 Google Sheet 同步。")
            return

        client = gspread.authorize(creds)
        try:
            spreadsheet = client.open(sheet_name)
        except gspread.SpreadsheetNotFound:
            spreadsheet = client.create(sheet_name)

        existing_worksheets = {ws.title: ws for ws in spreadsheet.worksheets()}

        for virus_name, df in virus_df_dict.items():
            ws_title = virus_name[:100]
            data_to_write = [df.columns.tolist()] + df.fillna("").astype(str).values.tolist() if not df.empty else []

            if data_to_write:
                if ws_title in existing_worksheets:
                    worksheet = existing_worksheets[ws_title]
                    worksheet.clear()
                    worksheet.update(values=data_to_write, range_name="A1")
                else:
                    worksheet = spreadsheet.add_worksheet(title=ws_title, rows=str(len(df)+10), cols="20")
                    worksheet.update(values=data_to_write, range_name="A1")

        logging.info(f"Google Sheet '{sheet_name}' 更新完畢！連結: {spreadsheet.url}")
    except Exception as e:
        logging.error(f"寫入 Google Sheet 發生錯誤: {e}")


# =============================================================================
# 核心任務執行 logic
# =============================================================================
def run_routine_search(num_sample_viruses: int = 5):
    # 1. 隨機選取病毒
    selected_viruses = random.sample(ANTIVIRALDB_RNA_VIRUS_LIST, k=min(num_sample_viruses, len(ANTIVIRALDB_RNA_VIRUS_LIST)))
    logging.info(f"=== 本次排程隨機抽樣病毒 ({len(selected_viruses)}種): {selected_viruses} ===")

    virus_to_cids: Dict[str, Set[int]] = {}
    virus_gene_info: Dict[str, Dict[str, str]] = {}
    all_unique_cids: Set[int] = set()

    # 2. 檢索 NCBI Gene 資訊與 PubChem CIDs
    for virus in selected_viruses:
        logging.info(f"--> [NCBI Gene] 查詢 '{virus}' 的 RdRP Domain 基因資訊...")
        gene_info = fetch_ncbi_gene_rdrp_info(virus)
        virus_gene_info[virus] = gene_info

        query_terms = [
            f"{virus} RdRP inhibitor",
            f"{virus} RNA-dependent RNA polymerase inhibitor"
        ]
        found_cids: Set[int] = set()
        for query in query_terms:
            cids = search_cids_via_esearch(query) or search_cids_via_pug_name(query)
            found_cids.update(cids)
            time.sleep(0.2)

        virus_to_cids[virus] = found_cids
        all_unique_cids.update(found_cids)
        logging.info(f"    找到 {len(found_cids)} 個抑制劑 CID | Gene Symbol: {gene_info['NCBI Gene Symbol']}")

    # 3. 批量獲取 PubChem 屬性與 SMILES 補全
    compounds_data = annotator.fetch_compound_properties_batch(list(all_unique_cids))
    pubchem_dict = {comp.get("CID"): comp for comp in compounds_data if comp.get("CID")}

    column_order = [
        "PubChem CID",
        "TargetViruses",
        "NCBI Gene ID",
        "NCBI Gene Symbol",
        "NCBI Gene Description",
        "Compound Pref Name",
        "Canonical SMILES by pubchem",
        "InChIKey",
        "Molecular Formula",
        "Molecular Weight"
    ]

    virus_df_dict: Dict[str, pd.DataFrame] = {}
    all_rows = []

    for virus in selected_viruses:
        cids = virus_to_cids.get(virus, set())
        gene_meta = virus_gene_info.get(virus, {})
        rows = []
        for cid in cids:
            comp = pubchem_dict.get(cid, {})
            comp_info = {
                "PubChem CID": cid,
                "TargetViruses": virus,
                "NCBI Gene ID": gene_meta.get("NCBI Gene ID", "N/A"),
                "NCBI Gene Symbol": gene_meta.get("NCBI Gene Symbol", "N/A"),
                "NCBI Gene Description": gene_meta.get("NCBI Gene Description", "N/A"),
                "Compound Pref Name": comp.get("Title", ""),
                "Canonical SMILES by pubchem": comp.get("Canonical SMILES by pubchem", ""),
                "InChIKey": comp.get("InChIKey", ""),
                "Molecular Formula": comp.get("MolecularFormula", ""),
                "Molecular Weight": comp.get("MolecularWeight", "")
            }
            rows.append(comp_info)
            all_rows.append(comp_info)

        df_v = pd.DataFrame(rows).reindex(columns=column_order) if rows else pd.DataFrame(columns=column_order)
        virus_df_dict[virus] = df_v

    df_combined = pd.DataFrame(all_rows).reindex(columns=column_order) if all_rows else pd.DataFrame(columns=column_order)

    # 存檔
    output_dir = "data"
    os.makedirs(output_dir, exist_ok=True)
    today_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_filename = os.path.join(output_dir, f"rdrp_inhibitors_random_{today_str}.csv")
    df_combined.to_csv(csv_filename, index=False, encoding="utf-8-sig")
    logging.info(f"本地 CSV 已儲存: {csv_filename}")

    # 上傳至 Google Sheets
    update_google_sheet_by_virus(virus_df_dict)


def main():
    INTERVAL_HOURS = 12  # 每 12 小時 (半天) 執行一次
    
    while True:
        logging.info(f"\n=================== 開始執行排程任務 ({datetime.now().strftime('%Y-%m-%d %H:%M:%S')}) ===================")
        try:
            # 每次隨機抽出 5 種病毒查詢
            run_routine_search(num_sample_viruses=5)
        except Exception as e:
            logging.error(f"執行任務時發生未預期錯誤: {e}")

        logging.info(f"=== 任務完成，進入休眠。下次執行時間為 {INTERVAL_HOURS} 小時後 ===")
        time.sleep(INTERVAL_HOURS * 3600)


if __name__ == "__main__":
    main()
