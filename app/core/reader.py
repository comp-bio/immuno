import re
from core.biogl import fasta_parse, translate, rev_comp


def is_dna(seq):
    return re.fullmatch(r'[ACGTNRYWSMKHBVDU]*', seq) is not None


def is_protein(seq):
    return re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYBXZ\*]*', seq) is not None


def src_parse(src):
    nuc, aac = (0, 0)
    sequences = []
    for name, s in fasta_parse(src, separator=""):
        seq = s.strip().upper().replace('-', '')
        sequences.append({'name': name, 'seq': seq})
        if is_dna(seq):
            nuc += 1
        elif is_protein(seq):
            aac += 1

    if aac > nuc:
        # Очистка от символов, не входящих в алфавит белка
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        regex = f"[^{re.escape(alphabet)}]"
        for item in sequences:
            item['seq'] = re.sub(regex, "*", item['seq'])
        return sequences
    else:
        return translate_it(sequences)


def translate_it(sequences):
    translated = []
    for item in sequences:
        s, name = (item['seq'], item['name'])
        for phase in range(0, 3):
            translated.append({
                'name': f"{name}___p{phase}-F",
                'seq': translate(s, phase=phase).replace('X', '*')
            })
        s = rev_comp(s)
        for phase in range(0, 3):
            translated.append({
                'name': f"{name}___p{phase}-R",
                'seq': translate(s, phase=phase).replace('X', '*')
            })
    return translated


def splitter(sequences, c=0.5):
    SIZES = [8, 9, 10, 11, 12, 13, 14]
    all_peps = {}
    for size in SIZES:
        for p in sequences:
            if 'chop' not in p:
                continue
            for i in range(0, len(p['seq']) - size):
                if p['chop'][i+size] < c:
                    continue
                one = p['seq'][i:i+size]
                if '*' in one:
                    continue
                all_peps[one] = True
    return [p for p in all_peps]
