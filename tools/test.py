import requests
import json

url = "https://neo-bolide.ndc.nasa.gov/service/event/public?limit=3"
r = requests.get(url)
print(r.status_code)
data = r.json()
print("Events found:", len(data.get('data', [])))
print(json.dumps(data.get('data', [])[0], indent=2)[:500])  # first event snippet