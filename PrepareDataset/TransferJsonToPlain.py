import json

def main():
    with open("./predicate_valuestype.json") as f:
        with open("./HopIndexPredicates.txt","w") as fout:
            s = f.read()
            a = json.loads(s)
            bindings = a['results']['bindings']
            for bind in bindings:
                fout.write("<"+bind['p']['value']+">\n")
main()