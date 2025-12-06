#!/usr/bin/env python3
"""Introspect gnomAD GraphQL schema."""

import requests
import json

gnomad_url = "https://gnomad.broadinstitute.org/api/"

# GraphQL introspection query for TranscriptConsequence type
introspection_query = """
query IntrospectionQuery {
    __type(name: "TranscriptConsequence") {
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

print("Introspecting gnomAD GraphQL schema...")
print("Looking for TranscriptConsequence fields...")
print()

try:
    response = requests.post(
        gnomad_url,
        json={"query": introspection_query},
        headers={"Content-Type": "application/json"},
        timeout=30
    )

    if response.status_code == 200:
        data = response.json()
        if "data" in data and "__type" in data["data"]:
            type_info = data["data"]["__type"]
            if type_info:
                print(f"Type: {type_info['name']}")
                print(f"\nAvailable fields:")
                print("-" * 60)
                for field in sorted(type_info["fields"], key=lambda x: x["name"]):
                    field_name = field["name"]
                    field_type = field["type"]["name"] or field["type"]["ofType"]["name"]
                    print(f"  {field_name:30s} : {field_type}")
            else:
                print("TranscriptConsequence type not found in schema")
        else:
            print(f"Unexpected response: {json.dumps(data, indent=2)[:500]}")
    else:
        print(f"Error: {response.status_code}")
        print(response.text[:500])

except Exception as e:
    print(f"Error: {e}")


# Also check Variant type
print("\n" + "="*60)
print("Looking for Variant fields...")
print("="*60)

variant_introspection = """
query IntrospectionQuery {
    __type(name: "Variant") {
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

try:
    response = requests.post(
        gnomad_url,
        json={"query": variant_introspection},
        headers={"Content-Type": "application/json"},
        timeout=30
    )

    if response.status_code == 200:
        data = response.json()
        if "data" in data and "__type" in data["data"]:
            type_info = data["data"]["__type"]
            if type_info:
                print(f"\nType: {type_info['name']}")
                print(f"\nAvailable fields (showing fields containing 'transcript'):")
                print("-" * 60)
                for field in sorted(type_info["fields"], key=lambda x: x["name"]):
                    field_name = field["name"]
                    if "transcript" in field_name.lower() or field_name in ["variant_id", "pos", "ref", "alt", "exome", "genome"]:
                        field_type = field["type"]["name"] or (field["type"].get("ofType") or {}).get("name", "complex")
                        print(f"  {field_name:30s} : {field_type}")

except Exception as e:
    print(f"Error: {e}")
