import uvicorn

from seqdelta.web import app


if __name__ == "__main__":
    uvicorn.run("seqdelta.web:app", host="127.0.0.1", port=8000, reload=False)
