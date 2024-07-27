import json

'''
Because datasets such as YAGO have multi-type problems, i.e. an entity
may have more than one type. But this problem may introduce conflict to
our assumption that each entity has only one type. 

So we must pre-process the dataset. For each entity having more than one
type with min number
'''

# record each type and it occur time.
type_num = {}
# the type each entity should have
result_type = {}


def load_type_num():
    # getAllType.json is the result of getAllType.sql
    with open("getAllType.json", encoding='utf-8') as f:
        t = json.load(f)
        result_list = t['results']['bindings']
        for tuple_result in result_list:
            type_name = tuple_result['o']['value']
            type_number = int(tuple_result['time']['value'])
            type_num[type_name] = type_number


def load_dup_num():
    with open("dupType.json", encoding='utf-8') as f:
        with open("deleteYAGO.sql", 'w', encoding='utf-8')  as delete_f:
            count = 0
            delete_f.write("DELETE DATA {\n")
            t = json.load(f)
            result_list = t['results']['bindings']
            for tuple_result in result_list:
                e_name = tuple_result['e']['value']
                new_name = tuple_result['t1']['value']
                if e_name not in result_type:
                    result_type[e_name] = new_name
                else:
                    old_type = result_type[e_name]
                    old_num = type_num[old_type]
                    new_num = type_num[new_name]
                    count = count + 1
                    # replace the old
                    if new_num < old_num:
                        result_type[e_name] = new_name
                        delete_f.write(
                            "<%s> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <%s> .\n" % (e_name, old_type))
                    # remains the same
                    else:
                        delete_f.write(
                            "<%s> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <%s>. \n" % (e_name, new_name))
            delete_f.write("}\n")
        print("output %d delete sentence" % count)


def look_dup_num():
    with open("./dupType.json", encoding='utf-8') as f:
        t = json.load(f)
        result_list = t['results']['bindings']
        counter = 0
        for tuple_result in result_list:
            e_name = tuple_result['e']['value']
            type_name = tuple_result['t1']['value']
            num = type_num[type_name]
            counter = counter + 1
            print("line %d" % counter, e_name, type_name, " type num:%d"%num)
            if counter == 30:
                break


def main():
    load_type_num()
    # look_dup_num()
    # return
    load_dup_num()


main()
