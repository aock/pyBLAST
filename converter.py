import argparse
import glob
import requests
from pprint import pprint
import re
import time
import json
import csv

PROGRAM = 'blastn&MEGABLAST=on'

# DATABASE = 'nt+nr'
DATABASE = 'nt'

datafolder = 'data'

base_uri = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'

# regex rid and time
pattern_rid = r'RID\s[=]\s(.{11})'
prog_rid = re.compile(pattern_rid)
pattern_estimated_time = r'RTOE\s[=]\s(\d*)'
prog_et = re.compile(pattern_estimated_time)

# regex status
pattern_status_waiting = r'\s+Status=WAITING'
pattern_status_failed = r'\s+Status=FAILED'
pattern_status_unknown = r'\s+Status=UNKNOWN'
pattern_status_ready = r'\s+Status=READY'

prog_status_waiting = re.compile(pattern_status_waiting)
prog_status_failed = re.compile(pattern_status_failed)
prog_status_unknown = re.compile(pattern_status_unknown)
prog_status_ready = re.compile(pattern_status_ready)

STATUS = {
    0: 'WAITING',
    1: 'FAILED',
    2: 'UNKNOWN',
    3: 'READY',
    4: 'NOMATCH'
}

# complete dna sequence
pattern_query_key = r'<QueryKey>(\d+)<\/QueryKey>'
pattern_web_env = r'<WebEnv>(\S+)<\/WebEnv>'


def get_rid( dna_seq ):


    content = 'CMD=Put&PROGRAM=%s&DATABASE=%s&QUERY=%s' % (PROGRAM, DATABASE, dna_seq)
    request = base_uri + '?' + content
    response = requests.get(request)

    rid_match = prog_rid.search(response.text)
    et_match = prog_et.search(response.text)

    rid = rid_match.group(1)
    print("rid: " + rid)

    et = et_match.group(1)
    print("estimated time: " + et)

    return rid, int(et)

def check_status( rid ):
    content = 'CMD=Get&FORMAT_OBJECT=SearchInfo&RID=%s' % rid
    request = base_uri + '?' + content
    response = requests.get(request)

    if prog_status_ready.search(response.text):
        print("ready to fetch!")
        return 3
    elif prog_status_waiting.search(response.text):
        print("waiting...")
        return 0
    elif prog_status_failed.search(response.text):
        print("failed!")
        return 1
    elif prog_status_unknown.search(response.text):
        print("search expired!")
        return 2
    else:
        print("No regex hits found?")
        return 4

def get_results( rid ):
    content = 'CMD=Get&FORMAT_TYPE=XML&RID=%s' % rid
    request = base_uri + '?' + content
    response = requests.get(request)
    return response.text

def name_dna(filename):

    name = ''
    dna = ''

    with open(filename, 'r') as in_file:
        for line in in_file:
            line = line.strip()

            if len(line) > 0 and line[0] == '>':
                if dna != '':
                    yield name, dna
                dna = ''
                name = line
            else:
                dna += line

    if dna != '':
        yield name,dna

def extract_best_dna(xml_string):
    # return name, accession
    # TODO
    return 'Duganella sp. c1 16S ribosomal RNA gene, partial sequence', 'JQ745646'

def get_complete_dna(hit_accession):
    # hit accession is code. e.g.: JQ745646
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=JQ745646&usehistory=y'

    request = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=%s&usehistory=y&retmode=json' % hit_accession
    response = requests.get(request)
    data = response.json()

    webenv = data['esearchresult']['webenv']
    query_key = data['esearchresult']['querykey']

    request = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&query_key=%s&rettype=fasta&retmode=text&WebEnv=%s' % (query_key, webenv)
    response = requests.get(request)

    fasta = response.text

    res = fasta.split('\n', 1)
    name = res[0]
    dna = res[1]

    return name, dna

def main():

    parser = argparse.ArgumentParser(description='Process Options.')
    parser.add_argument('-i','--input', type=str, help='input filename', required=True)
    parser.add_argument('-o','--output', type=str, help='output filename', required=True)
    args = parser.parse_args()


    retry_count = 20

    table = []

    for name, dna in name_dna(args.input):
        print(name)
        print(dna)

        # first get request id with dna sequence
        rid, et = get_rid(dna)

        # wait for search to complete
        time.sleep(et)

        counter = 0
        status = 0
        while status != 3:
            time.sleep(5)
            print("fetch result try: " + str(counter))
            status = check_status(rid)
            counter += 1

            if counter >= retry_count:
                print("retry")
                rid, et = get_rid(dna)
                time.sleep(et)
                counter = 0

        res = get_results(rid)

        # extract best dna -> hit_accession
        a_name, hit_accession = extract_best_dna(res)
        c_name, c_dna = get_complete_dna(hit_accession)
        table.append([name, dna, c_name, c_dna])

    with open("output.csv", "wb") as f:
        writer = csv.writer(f)
        writer.writerows(table)



if __name__ == "__main__":
    main()