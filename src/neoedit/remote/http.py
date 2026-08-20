"""Minimal HTTP helper on urllib: user agent, timeouts, polite retries on 429/5xx."""
from __future__ import annotations

import json
import os
import re
import socket
import time
import urllib.error
import urllib.parse
import urllib.request

from .. import __version__

USER_AGENT = f"NeoEdit/{__version__} (+https://github.com/evozoa/NeoEdit)"


class RemoteError(Exception):
    """A request failed in a way worth showing to the user (network, HTTP error, server message)."""


def _host(url: str) -> str:
    return urllib.parse.urlsplit(url).netloc or url


def _server_message(body: str) -> str:
    """Pull a human-readable message out of an error body (Ensembl JSON, NCBI text, HTML)."""
    b = body.strip()
    if b.startswith("{"):
        try:
            d = json.loads(b)
            if isinstance(d, dict) and d.get("error"):
                return str(d["error"])
        except ValueError:
            pass
    if "<html" in b[:400].lower():
        m = re.search(r"<title>(.*?)</title>", b, re.S | re.I)
        return (m.group(1).strip() if m else "HTML error page")
    return b[:300]


def http_get(url: str, *, headers: dict | None = None, data: bytes | None = None,
             timeout: float = 60.0, retries: int = 3) -> bytes:
    """GET (or POST when `data` is given) and return the body. Retries transient errors.

    Raises RemoteError with a readable message on failure."""
    hdrs = {"User-Agent": USER_AGENT, **(headers or {})}
    last_err = None
    for attempt in range(retries + 1):
        req = urllib.request.Request(url, data=data, headers=hdrs)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            body = ""
            try:
                body = e.read().decode("utf-8", "replace")
            except Exception:
                pass
            if e.code in (429, 500, 502, 503, 504) and attempt < retries:
                try:
                    wait = float(e.headers.get("Retry-After", "") or (1.0 + attempt))
                except ValueError:
                    wait = 1.0 + attempt
                time.sleep(min(wait, 30.0))
                last_err = e
                continue
            msg = _server_message(body)
            raise RemoteError(f"{_host(url)} answered HTTP {e.code}" + (f": {msg}" if msg else "")) from None
        except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError) as e:
            reason = getattr(e, "reason", e)
            if attempt < retries and not isinstance(reason, (socket.gaierror,)):
                time.sleep(1.0 + attempt)
                last_err = e
                continue
            raise RemoteError(f"Cannot reach {_host(url)}: {reason}") from None
    raise RemoteError(f"Giving up on {_host(url)}: {last_err}")


def safe_filename(stem: str, maxlen: int = 80) -> str:
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", stem).strip("._")
    return (s or "download")[:maxlen]


def write_download(text: str, out_dir: str, stem: str, ext: str) -> str:
    """Write `text` to out_dir/<safe stem><ext>, numbering the name if a *different* file is there."""
    os.makedirs(out_dir, exist_ok=True)
    base = safe_filename(stem)
    path = os.path.join(out_dir, base + ext)
    n = 1
    while os.path.exists(path):
        try:
            with open(path, "r", errors="replace") as fh:
                if fh.read() == text:
                    return path                       # identical content already there: reuse it
        except OSError:
            pass
        n += 1
        path = os.path.join(out_dir, f"{base}_{n}{ext}")
    with open(path, "w") as fh:
        fh.write(text)
    return path
