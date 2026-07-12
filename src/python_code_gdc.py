def get_gdc_clinical_supplements(
    case_ids: pd.Series,
    batch_size=200,
):
    endpoint = "https://api.gdc.cancer.gov/files"

    case_ids = np.unique(case_ids)

    columns = [
        "case_id",
        "barcode_case",
        "project_id",
        "file_id",
        "file_name",
        "data_category",
        "data_type",
        "data_format",
        "access",
    ]

    records = []
    icount=-1
    batchcount=-1

    for batch_number, start in enumerate(
        range(0, len(case_ids), batch_size),
        start=1,
    ):
        batch = case_ids[start:start + batch_size]

        batchcount += 1
        print(f"b{batchcount}",end=' ')

        if isinstance(batch, np.ndarray):
            batch = batch.tolist()    

        print(";".join(batch))   

        filters = {
            "op": "and",
            "content": [
                {
                    "op": "in",
                    "content": {
                        "field": "cases.case_id",
                        "value": batch,
                    },
                },
            ],
        }

        params = {
            "filters": json.dumps(filters),
            "fields": ",".join(
                [
                    "file_id",
                    "file_name",
                    "data_category",
                    "data_type",
                    "data_format",
                    "access",
                    "cases.case_id",
                    "cases.submitter_id",
                    "cases.project.project_id",
                ]
            ),
            "format": "JSON",
            "size": 10000,
        }

        try:
            response = requests.get(
                endpoint,
                params=params,
                timeout=120,
            )
            response.raise_for_status()
        except requests.RequestException as exc:
            print(
                f"\nBatch {batch_number} failed: {exc}"
            )
            continue

        payload = response.json()

        hits = payload.get("data", {}).get("hits", [])
        pagination = payload.get("data", {}).get(
            "pagination",
            {},
        )

        total = pagination.get("total", len(hits))

        print(
            f"Batch {batch_number}: "
            f"{len(batch)} cases, "
            f"{total} clinical files"
        )

        for hit in hits:
            icount += 1

            if icount % 10 == 0:
                print(".", end="", flush=True)

            associated_cases = hit.get("cases") or []

            if not associated_cases:
                records.append(
                    {
                        "case_id": None,
                        "barcode_case": None,
                        "project_id": None,
                        "file_id": hit.get("file_id"),
                        "file_name": hit.get("file_name"),
                        "data_category": hit.get(
                            "data_category"
                        ),
                        "data_type": hit.get("data_type"),
                        "data_format": hit.get(
                            "data_format"
                        ),
                        "access": hit.get("access"),
                    }
                )
                continue

            for case in associated_cases:
                project = case.get("project") or {}

                records.append(
                    {
                        "case_id": case.get("case_id"),
                        "barcode_case": case.get(
                            "submitter_id"
                        ),
                        "project_id": project.get(
                            "project_id"
                        ),
                        "file_id": hit.get("file_id"),
                        "file_name": hit.get("file_name"),
                        "data_category": hit.get(
                            "data_category"
                        ),
                        "data_type": hit.get("data_type"),
                        "data_format": hit.get(
                            "data_format"
                        ),
                        "access": hit.get("access"),
                    }
                )

    print("\n--------------- end ------------")

    if records == []:
        df = pd.DataFrame()
    else:
        try:
            df = pd.DataFrame(records)
        except Exception as e:
            print(f"Error: {e}")
            df = pd.DataFrame()

    return df


#===========================

def verify_gdc_cases(case_ids):
    endpoint = "https://api.gdc.cancer.gov/cases"

    case_ids = (
        pd.Series(case_ids, dtype="object")
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    filters = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [gdc.gdc_project_id],
                },
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.primary_site",
                    "value": ["Breast"],
                },
            },
        ],
    }

    params = {
        "filters": json.dumps(filters),
        "fields": ",".join(
            [
                "case_id",
                "submitter_id",
                "project.project_id",
                "primary_site",
                "disease_type",
            ]
        ),
        "format": "JSON",
        "size": len(case_ids),
    }

    response = requests.get(
        endpoint,
        params=params,
        timeout=120,
    )
    response.raise_for_status()

    hits = response.json()["data"]["hits"]

    records = []

    for case in hits:
        project = case.get("project") or {}

        records.append(
            {
                "case_id": case.get("case_id"),
                "barcode_case": case.get("submitter_id"),
                "project_id": project.get("project_id"),
                "primary_site": case.get("primary_site"),
                "disease_type": case.get("disease_type"),
            }
        )

    return pd.DataFrame(records)

df_verify = verify_gdc_cases(df_cases["case_id"].head(50))

print(df_verify.shape)
print(df_verify["project_id"].value_counts(dropna=False))
display(df_verify.head())


#===========================




filters = {
    "op": "and",
    "content": [
        {
            "op": "in",
            "content": {
                "field": "cases.case_id",
                "value": df_cases["case_id"].tolist(),
            },
        },
        {
            "op": "in",
            "content": {
                "field": "data_category",
                "value": ["Transcriptome Profiling"],
            },
        },
        {
            "op": "in",
            "content": {
                "field": "data_type",
                "value": ["Gene Expression Quantification"],
            },
        },
    ],
}

params = {
    "filters": json.dumps(filters),
    "fields": ",".join(
        [
            "file_id",
            "file_name",
            "data_category",
            "data_type",
            "data_format",
            "workflow_type",
            "analysis.workflow_type",
            "cases.case_id",
            "cases.submitter_id",
            "cases.samples.sample_id",
            "cases.samples.submitter_id",
            "cases.samples.sample_type",
        ]
    ),
    "format": "JSON",
    "size": 10000,
}

response = requests.get(
    "https://api.gdc.cancer.gov/files",
    params=params,
    timeout=120,
)
response.raise_for_status()

hits = response.json()["data"]["hits"]

len(hits)

expression_records = []

for hit in hits:
    analysis = hit.get("analysis") or {}

    for case in hit.get("cases") or []:
        expression_records.append(
            {
                "case_id": case.get("case_id"),
                "barcode_case": case.get("submitter_id"),
                "file_id": hit.get("file_id"),
                "file_name": hit.get("file_name"),
                "data_format": hit.get("data_format"),
                "workflow_type": (
                    hit.get("workflow_type")
                    or analysis.get("workflow_type")
                ),
            }
        )

df_expression_files = (
    pd.DataFrame(expression_records)
    .drop_duplicates()
    .reset_index(drop=True)
)

df_expression_files


