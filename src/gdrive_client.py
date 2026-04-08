#!/usr/bin/python
#!python
# -*- coding: utf-8 -*-
# Created on 2026/04/07
# Updated  on 2026/04/07
# @author: Flavio Lichtenstein
# @local: Home sweet home

import io
import os
import json
from pathlib import Path
from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload

SCOPES = ["https://www.googleapis.com/auth/drive.readonly"]

def get_drive_service():
    raw = os.environ["GOOGLE_SERVICE_ACCOUNT_JSON"]
    info = json.loads(raw)
    creds = service_account.Credentials.from_service_account_info(
        info,
        scopes=SCOPES,
    )
    return build("drive", "v3", credentials=creds)

def list_folder_files(folder_id: str):
    service = get_drive_service()
    query = f"'{folder_id}' in parents and trashed = false"
    resp = service.files().list(
        q=query,
        fields="files(id,name,mimeType,modifiedTime,size)",
        includeItemsFromAllDrives=True,
        supportsAllDrives=True,
    ).execute()
    return resp.get("files", [])

def download_file(file_id: str, dest_path: str):
    service = get_drive_service()
    request = service.files().get_media(fileId=file_id)

    Path(dest_path).parent.mkdir(parents=True, exist_ok=True)
    with open(dest_path, "wb") as fh:
        downloader = MediaIoBaseDownload(fh, request)
        done = False
        while not done:
            _, done = downloader.next_chunk()