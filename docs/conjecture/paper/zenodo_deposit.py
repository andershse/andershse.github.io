"""Preview, create, or publish the prepared software deposit on Zenodo.

Uses only the Python standard library. A token is read from a hidden terminal
prompt or ZENODO_ACCESS_TOKEN, sent only in HTTPS headers, and never saved.
Default action is an offline preview. The explicit 'publish' action uploads
and publishes the archive with the metadata/license shown in that preview.
"""

import argparse
from datetime import date
import getpass
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time
from urllib.error import HTTPError
from urllib.parse import quote, urlparse
from urllib.request import HTTPRedirectHandler, Request, build_opener
import zipfile

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "identifiability-verification.zip"
STATE = HERE / "zenodo_record.json"
BASE = "https://zenodo.org"
REQUIRED = {
    "README.txt", "LICENSE", "CITATION.cff", "requirements.txt",
    "verify_counterexample.py", "verify_spectral_counterexample.py",
    "verification_results.json", "spectral_verification_results.json",
    "MANIFEST.json",
}


def verify_archive(path):
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        if len(names) != len(REQUIRED) or set(names) != REQUIRED:
            raise ValueError("The archive contains missing, duplicate, or unexpected files.")
        manifest = json.loads(archive.read("MANIFEST.json"))
        if set(manifest) != REQUIRED - {"MANIFEST.json"}:
            raise ValueError("The archive manifest has unexpected entries.")
        for name, expected in manifest.items():
            if hashlib.sha256(archive.read(name)).hexdigest() != expected:
                raise ValueError("Archive checksum mismatch for " + name)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_url(url):
    parsed = urlparse(url)
    if (parsed.scheme, parsed.netloc) != ("https", "zenodo.org"):
        raise ValueError("Refusing to send credentials outside https://zenodo.org.")


class ZenodoRedirects(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        check_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def api(token, method, url, payload=None, file_data=None):
    check_url(url)
    headers = {"Authorization": "Bearer " + token, "Accept": "application/json"}
    body = file_data
    if payload is not None:
        body = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    elif file_data is not None:
        headers["Content-Type"] = "application/zip"
    request = Request(url, data=body, headers=headers, method=method)
    with build_opener(ZenodoRedirects()).open(request, timeout=60) as response:
        return json.load(response)


def save_state(record, checksum):
    metadata = record.get("metadata", {})
    info = {
        "deposition_id": record["id"],
        "archive_sha256": checksum,
        "draft_url": record.get("links", {}).get("html"),
        "doi": metadata.get("doi") or
               metadata.get("prereserve_doi", {}).get("doi"),
        "submitted": bool(record.get("submitted")),
    }
    STATE.write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    return info


def write_citation(record):
    doi = record.get("doi") or record.get("metadata", {}).get("doi")
    if not doi or not re.fullmatch(r"10\.5281/zenodo\.\d+", doi):
        raise ValueError("A confirmed production Zenodo DOI was not returned.")
    (HERE / "verification_doi.tex").write_text(
        "\\newcommand{\\verificationdoi}{" + doi + "}\n", encoding="utf-8")
    print("Public software record: https://doi.org/" + doi)
    print("The manuscript DOI include has been written. Recompile the paper to add it.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("preview", "draft", "publish"),
                        default="preview", nargs="?")
    args = parser.parse_args()
    checksum = verify_archive(ARCHIVE)
    payload = json.loads((HERE / "zenodo_metadata.json").read_text(encoding="utf-8"))
    payload["metadata"]["publication_date"] = date.today().isoformat()
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    print("Archive:", ARCHIVE)
    print("SHA-256:", checksum)
    if args.action == "preview":
        print("Offline preview only: no account was accessed and nothing was uploaded.")
        return

    token = os.environ.get("ZENODO_ACCESS_TOKEN")
    if not token:
        if not sys.stdin.isatty():
            raise ValueError("Run in an interactive terminal to enter the token securely.")
        token = getpass.getpass("Zenodo token (hidden; not saved): ").strip()
    if not token:
        raise ValueError("A Zenodo token is required.")

    try:
        if STATE.exists():
            prior = json.loads(STATE.read_text(encoding="utf-8"))
            identifier = int(prior["deposition_id"])
            record = api(token, "GET", BASE + f"/api/deposit/depositions/{identifier}")
            if record.get("submitted"):
                if checksum != prior["archive_sha256"]:
                    raise ValueError("The published archive differs. Create a new version in Zenodo.")
            else:
                print("Resuming draft:", identifier)
        else:
            record = api(token, "POST", BASE + "/api/deposit/depositions", payload={})
            save_state(record, checksum)

        endpoint = BASE + f"/api/deposit/depositions/{int(record['id'])}"
        if not record.get("submitted"):
            record = api(token, "PUT", endpoint, payload=payload)
            save_state(record, checksum)
            bucket = record["links"]["bucket"].rstrip("/")
            contents = ARCHIVE.read_bytes()
            upload = api(token, "PUT", bucket + "/" + quote(ARCHIVE.name),
                         file_data=contents)
            # Zenodo documents MD5 for the file-upload checksum.
            if upload.get("checksum", "").removeprefix("md5:") != hashlib.md5(contents).hexdigest():
                raise ValueError("Zenodo's uploaded-file checksum does not match.")
            record = api(token, "GET", endpoint)
            info = save_state(record, checksum)
            if args.action == "draft":
                print("Draft saved:", info["draft_url"])
                print("Reserved DOI:", info["doi"], "(not yet public)")
                return
            record = api(token, "POST", endpoint + "/actions/publish")
            save_state(record, checksum)

        public_endpoint = BASE + f"/api/records/{int(record.get('record_id', record['id']))}"
        for attempt in range(5):
            try:
                public = api(token, "GET", public_endpoint)
                write_citation(public)
                return
            except HTTPError as error:
                if error.code != 404:
                    raise
                if attempt < 4:
                    time.sleep(2)
        print("Publication was requested; the public record is still being processed.")
        print("Rerun this command to retrieve the DOI after Zenodo finishes.")
    except HTTPError as error:
        # Avoid printing request headers or the token, including on API errors.
        detail = error.read().decode("utf-8", errors="replace").replace(token, "[redacted]")
        raise RuntimeError(f"Zenodo returned HTTP {error.code}: {detail[:1200]}") from None


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError) as error:
        print("Deposit stopped:", str(error), file=sys.stderr)
        sys.exit(1)
