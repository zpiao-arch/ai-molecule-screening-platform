#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PORT="${PORT:-8765}"
HEALTH_URL="http://localhost:${PORT}/api/health"
SYSTEM_HEALTH_URL="http://localhost:${PORT}/api/system/health"
APP_URL="http://localhost:${PORT}/"

cd "$ROOT_DIR"

echo "Default local URL: http://localhost:8765/"
echo "Web app: ${APP_URL}"
echo "Health check: ${HEALTH_URL}"
echo "Product health: ${SYSTEM_HEALTH_URL}"
echo "Local health summary: python3 health_check.py"

if command -v lsof >/dev/null 2>&1; then
  EXISTING_PID="$(lsof -tiTCP:${PORT} -sTCP:LISTEN 2>/dev/null | head -n 1 || true)"
  if [ -n "${EXISTING_PID}" ]; then
    echo "Port ${PORT} is already in use by PID ${EXISTING_PID}."
    echo "Checking existing service..."
    EXISTING_READY=""
    for _ in $(seq 1 8); do
      if curl -sS "${HEALTH_URL}"; then
        EXISTING_READY="1"
        break
      fi
      sleep 1
    done
    echo
    if [ -z "${EXISTING_READY}" ]; then
      echo "Existing service is not responding to ${HEALTH_URL}."
      echo "Use a different port: PORT=8766 ./start_web.sh"
      echo "Or stop PID ${EXISTING_PID} if it is stale, then retry ./start_web.sh."
      exit 1
    fi
    echo "Existing service became ready."
    echo "Checking Product Ops health..."
    curl -sS "${SYSTEM_HEALTH_URL}" || true
    echo
    echo "If this is the project server, open ${APP_URL}."
    exit 0
  fi
fi

echo "Starting FastAPI server with webapp/server.py ..."
python3 webapp/server.py &
SERVER_PID="$!"

for _ in $(seq 1 30); do
  if curl -sS "${HEALTH_URL}" >/dev/null 2>&1; then
    echo "Server is ready. PID: ${SERVER_PID}"
    curl -sS "${HEALTH_URL}"
    echo
    curl -sS "${SYSTEM_HEALTH_URL}" || true
    echo
    echo "Open ${APP_URL} and enter 产品指挥台 for the Stage 8 demo."
    echo "For a readable local summary, run: python3 health_check.py"
    wait "${SERVER_PID}"
    exit 0
  fi
  sleep 1
done

echo "Server did not pass /api/health within 30 seconds."
echo "Check the terminal output above, then retry: ./start_web.sh"
exit 1
