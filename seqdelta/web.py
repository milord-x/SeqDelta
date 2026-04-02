from __future__ import annotations

from pathlib import Path
from uuid import uuid4

from fastapi import FastAPI, File, HTTPException, Request, UploadFile
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from .mutation import analyze_records
from .parser import parse_fasta_text
from .report import (
    PROJECT_ROOT,
    RESOLVED_STATIC_DIR,
    RESOLVED_TEMPLATES_DIR,
    build_report_context,
    write_csv_report,
    write_html_report,
    write_json_report,
)

GENERATED_DIR = PROJECT_ROOT / "generated_reports"
GENERATED_DIR.mkdir(parents=True, exist_ok=True)

templates = Jinja2Templates(directory=str(RESOLVED_TEMPLATES_DIR))
app = FastAPI(title="SeqDelta", version="0.1.0")
app.mount("/static", StaticFiles(directory=str(RESOLVED_STATIC_DIR)), name="static")
app.mount("/generated_reports", StaticFiles(directory=str(GENERATED_DIR)), name="generated_reports")


@app.get("/", response_class=HTMLResponse)
async def landing_page(request: Request):
    return templates.TemplateResponse(
        "index.html",
        {
            "request": request,
            "report_title": "SeqDelta",
        },
    )


@app.post("/analyze", response_class=HTMLResponse)
async def analyze_uploads(
    request: Request,
    reference_fasta: UploadFile = File(...),
    mutant_fasta: UploadFile = File(...),
):
    try:
        reference_text = (await reference_fasta.read()).decode("utf-8")
        mutant_text = (await mutant_fasta.read()).decode("utf-8")
        reference_record = parse_fasta_text(reference_text, default_identifier="reference", source_name=reference_fasta.filename)
        mutant_record = parse_fasta_text(mutant_text, default_identifier="mutant", source_name=mutant_fasta.filename)
    except UnicodeDecodeError as exc:
        raise HTTPException(status_code=400, detail="Uploaded FASTA files must be UTF-8 text.") from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    result = analyze_records(reference_record, mutant_record)

    analysis_id = uuid4().hex
    html_name = f"{analysis_id}.html"
    json_name = f"{analysis_id}.json"
    csv_name = f"{analysis_id}.csv"
    write_html_report(result, GENERATED_DIR / html_name)
    write_json_report(result, GENERATED_DIR / json_name)
    write_csv_report(result, GENERATED_DIR / csv_name)

    context = build_report_context(
        result,
        standalone_report=False,
        report_title="SeqDelta Results",
        include_plotly_bundle=True,
        download_links={
            "html": f"/generated_reports/{html_name}",
            "json": f"/generated_reports/{json_name}",
            "csv": f"/generated_reports/{csv_name}",
        },
    )
    context["request"] = request
    return templates.TemplateResponse("results.html", context)


@app.get("/healthz")
async def healthcheck():
    return {"status": "ok"}


@app.get("/favicon.ico")
async def favicon():
    return RedirectResponse(url="/static/favicon.svg")
