import json

'''
total triple in origin dataset: 46093317
numbers in DupType.sql
    pair num: 16352423
    3243477 distinct entity
    delete triples number: 13,108,946
32984371 triples are left after deleting
successfully delete!
'''

result_type = {}
def load_dup_num():
    with open("dupType.json", encoding='utf-8') as f:
        t = json.load(f)
        result_list = t['results']['bindings']
        print("pair num:",len(result_list))

        s = set([tuple_result['e']['value'] for tuple_result in result_list])
        print("%d distinct entity"%len(s))

def main():
    load_dup_num()


main()
