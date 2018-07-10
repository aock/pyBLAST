import argparse
import glob
import requests
from pprint import pprint
import re
import time
import json
import csv
import xml.etree.cElementTree as ET
import os
import json

PROGRAM = 'blastn&MEGABLAST=on'

# high to low
sort_descending = True

# avoid words from config
# avoiding_words = ['[Cc]omplete genome',
#     '[Uu]ncultured',
#     '[Uu]nknown']

# score words from config
# score_words = {
#     'complete genome' : 0.1,
#     '16S.*RNA': 1.5,
#     '[Uu]ncultured': 0.0
# }


# DATABASE = 'nt+nr'
# DATABASE = 'nt'
# database from config

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


def get_rid( dna_seq, config ):

    DATABASE = config['database']

    content = 'CMD=Put&PROGRAM=%s&DATABASE=%s&QUERY=%s' % (PROGRAM, DATABASE, dna_seq)
    request = base_uri + '?' + content
    response = requests.get(request)

    rid_match = prog_rid.search(response.text)
    et_match = prog_et.search(response.text)

    rid = rid_match.group(1)
    et = et_match.group(1)

    return rid, int(et)

def check_status( rid ):
    content = 'CMD=Get&FORMAT_OBJECT=SearchInfo&RID=%s' % rid
    request = base_uri + '?' + content
    response = requests.get(request)

    if prog_status_ready.search(response.text):
        return 3
    elif prog_status_waiting.search(response.text):
        return 0
    elif prog_status_failed.search(response.text):
        return 1
    elif prog_status_unknown.search(response.text):
        return 2
    else:
        return 4

def get_results( rid ):
    content = 'CMD=Get&FORMAT_TYPE=XML&RID=%s' % rid
    request = base_uri + '?' + content
    response = requests.get(request)
    return response.text

def name_dna(filename):

    name = ''
    dna = ''
    dna_orig = ''

    with open(filename, 'r') as in_file:
        for line in in_file:
            line = line.strip()

            if len(line) > 0 and line[0] == '>':
                if dna != '':
                    yield name, dna, dna_orig.rstrip('\n')
                dna = ''
                dna_orig = ''
                name = line
            else:
                dna += line
                dna_orig += line + '\n'

    if dna != '':
        yield name,dna,dna_orig.rstrip('\n')

def extract_best_dna(xml_string, config):
    # return name, accession

    sort_kw = config['default_score']
    avoiding_words = config['avoid_words']
    score_words = config['score_words']

    hit_dicts = []

    root = ET.fromstring(xml_string)
    hits = root.find('BlastOutput_iterations').find('Iteration').find('Iteration_hits').findall('Hit')

    # generate list of dictionaries containing several important information about results
    for hit in hits:
        hit_dict = {}
        hit_dict['accession'] = hit.find('Hit_accession').text
        hit_dict['name'] = hit.find('Hit_def').text

        # AVOID
        skip = False
        for avoid_word in avoiding_words:
            m = re.search(avoid_word, hit_dict['name'])
            if m is not None:
                skip = True

        if skip:
            continue

        hit_dict['id'] = hit.find('Hit_id').text

        scores = hit.find('Hit_hsps').find('Hsp')

        hit_dict['bit-score'] = float(scores.find('Hsp_bit-score').text)
        hit_dict['score'] = int(scores.find('Hsp_score').text)
        hit_dict['identity'] = int(scores.find('Hsp_identity').text)
        hit_dict['hit-from'] = int(scores.find('Hsp_hit-from').text)
        hit_dict['hit-to'] = int(scores.find('Hsp_hit-to').text)

        # SCORE
        hit_dict['own_score'] = hit_dict[sort_kw]

        for k,v in score_words.items():
            m = re.search(k, hit_dict['name'])
            if m is not None:
                hit_dict['own_score'] *= v

        hit_dicts.append(hit_dict)

    # sorting list
    hit_dicts = sorted(hit_dicts, key=lambda k: k['own_score'], reverse = sort_descending)

    return hit_dicts[0]

def get_complete_dna(hit_accession):
    # hit accession is code. e.g.: JQ745646

    name = 'bla'
    dna = 'bla'

    while name[0] != '>':

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
        dna = res[1].strip('\n')


    return name, dna

def test_xml_parser(config):
    xml_str = ''
    with open('test.xml','r') as f:
        xml_str = f.read()

    name, accession = extract_best_dna(xml_str, config)

    print(name)
    print(accession)

def save_response(res):
    with open('test.xml','w') as f:
        f.write(res)

def read_response():
    with open('test.xml','r') as f:
        return f.read()

def main():

    parser = argparse.ArgumentParser(description='Process Options.')
    parser.add_argument('-i','--input', type=str, help='input filename', required=True)
    parser.add_argument('-o','--output', type=str, help='output filename', required=True)
    parser.add_argument('-c','--config', type=str, help='config file (default: settings.json)', required=False, default='settings.json')
    args = parser.parse_args()

    config = {}
    with open(args.config,'r') as f:
        config = json.load(f)

    pprint(config)


    retry_count = 40

    column_names = [
            'name',
            'dna_orig',
            'c_name',
            'c_dna',
            'accession',
            'own_score',
            'bit-score',
            'score',
            'identity',
            'hit-from',
            'hit-to']   
    table = []


    input_data = [(name,dna,dna_orig) for name, dna, dna_orig in name_dna(args.input) ]



    for i,(name, dna, dna_orig) in enumerate(input_data):

        print(str(i) + '/' + str(len(input_data)) + ': ' + name)

        while True:

            # first get request id with dna sequence
            rid, et = get_rid(dna, config)
            print(et)

            # wait for search to complete
            time.sleep(et)

            counter = 0
            status = 0
            while status != 3:
                time.sleep(5)
                status = check_status(rid)
                counter += 1

                if counter >= retry_count:
                    rid, et = get_rid(dna, config)
                    print(et)
                    time.sleep(et)
                    counter = 0

            res = get_results(rid)

            # only debug. comment this
            # save xml response
            # save_response(res)
            # print('saved')
            # exit()
            # read xml response
            # res = read_response()

            # extract best dna -> hit_accession
            try:
                hit = extract_best_dna(res, config)
                c_name, c_dna = get_complete_dna( hit['accession'] )

                table.append([name,
                            dna_orig,
                            c_name,
                            c_dna,
                            hit['accession'],
                            hit['own_score'],
                            hit['bit-score'],
                            hit['score'],
                            hit['identity'],
                            hit['hit-from'],
                            hit['hit-to']])

                break
            except:
                with open('error/error_'+name+'.xml','w') as f:
                    f.write(res)
                print("error parsing xml")


    out_folder = args.output.strip('/')

    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    
    # create csv with containing all information
    with open(args.output + '/conv.csv', 'w') as f:
        writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(table)

    # creating weird fasta format file
    with open(out_folder + '/fasta.txt', 'w') as f:
        for i,row in enumerate(table):
            for j in range(4):
                f.write(row[j])
                f.write('\n')
                f.write('\n')

    



if __name__ == "__main__":
    main()
