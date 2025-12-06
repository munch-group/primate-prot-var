#!/usr/bin/env python3
"""Debug gnomAD GraphQL query - introspect schema."""

import requests
import json

gnomad_url = "https://gnomad.broadinstitute.org/api/"

# Introspect Fafmax
query = """
query IntrospectionQuery {
    __type(name: "Fafmax") {
        name
        fields {
            name
            type {
                name
                kind
                ofType {
                    name
                    kind
                }
            }
        }
    }
}
"""

response = requests.post(
    gnomad_url,
    json={"query": query},
    headers={"Content-Type": "application/json"},
    timeout=30
)

print("=== Fafmax fields ===")
data = response.json()
if data.get("data", {}).get("__type"):
    for field in sorted(data["data"]["__type"]["fields"], key=lambda x: x["name"]):
        field_type = field["type"]["name"] or (field["type"].get("ofType") or {}).get("name", "complex")
        print(f"  {field['name']:30s} : {field_type}")

# Introspect VariantFilteringAlleleFrequency
query2 = """
query IntrospectionQuery {
    __type(name: "VariantFilteringAlleleFrequency") {
        name
        fields {
            name
            type {
                name
                kind
                ofType {
                    name
                    kind
                }
            }
        }
    }
}
"""

response2 = requests.post(
    gnomad_url,
    json={"query": query2},
    headers={"Content-Type": "application/json"},
    timeout=30
)

print("\n=== VariantFilteringAlleleFrequency fields ===")
data2 = response2.json()
if data2.get("data", {}).get("__type"):
    for field in sorted(data2["data"]["__type"]["fields"], key=lambda x: x["name"]):
        field_type = field["type"]["name"] or (field["type"].get("ofType") or {}).get("name", "complex")
        print(f"  {field['name']:30s} : {field_type}")
