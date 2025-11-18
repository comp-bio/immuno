
import json
import glob
import os
import csv
from config import log


def reader(src):
    objects = {}
    with open(src, 'r') as file:
        table = list(csv.reader(file, delimiter="\t"))
        header, data = (table[0], table[1:])
        for line in data:
            if not line:
                continue
            obj = dict(zip(header, line))
            if obj['Sequence'] not in objects:
                objects[obj['Sequence']] = []
            objects[obj['Sequence']].append(obj)
    return objects


def search(db, items):
    results = {}
    for db_name, db_data in db:
        results[db_name] = {}
        for prot in items:
            if prot not in db_data:
                continue
            results[db_name][prot] = db_data[prot]
    return results


def IEAtlas(params):
    for src in glob.glob(f"./db/IEAtlas/Epitopes*.txt"):
        yield (os.path.basename(src).replace('.txt', ''), reader(src))


def IEDB(params):
    for src in glob.glob(f"./db/epitope/{params.organism}.json"):
        data = {}
        with open(src) as js:
            content = json.load(js)['Data']
            for item in content:
                name = item['Epitope - Name']
                if name not in data:
                    data[name] = []
                data[name].append(item)
        yield (os.path.basename(src).replace('.json', ''), data)


def export_all(info):
    message = ""
    for file in info:
        if len(info[file]) > 0:
            message += ("-" * 80) + "\n"
            message += f"{file}:\n"
            message += ("-" * 80) + "\n"
        for pep in info[file]:
            for block in info[file][pep]:
                for k in block:
                    if block[k]:
                        message += f'{k}: {block[k]}\n'
            message += "\n"
        message += "\n"
    return message


def export_file(params, name, info):
    total = sum([len(v) for v in info.values()])
    if total == 0:
        log(f"– {name} ({total})")
        return

    message = export_all(info)
    out = f"{params.dir}/extdb-{name}.txt"
    with open(out, 'w+') as res:
        res.write(message)
    log(f"– {name} ({total}). File: {out}")
