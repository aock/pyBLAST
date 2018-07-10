# pyBLAST

Automatic matching of multiple DNA sequences with DNAs from BLAST databases via Web-REST interface. Additional functionalities to improve results:

<ul>
  <li> Define keywords, that are not allowed in resulting DNA sequence name </li>
  <li> Define keywords occurence with multiplicators </li>
</ul>

## Usage

```console
user@pc:~$ python converter.py --input data/example.txt -o dna1
```

input file format (FASTA) for "data/example.txt":

```
>unknown_name_1
CAGTCGAACGGCCA...
CAGTCGAACGGCCA...
...
>unknown_name_2
CAGTCGAACGGCCA...
CAGTCGAACGGCCA...
...
```

software generates folder named "dna1", containing "conv.csv" and "fasta.txt"

"fasta.txt" is a FASTA formatted file:

```
>unknown_name_1

CAGTCGAACGGCCA...
CAGTCGAACGGCCA...
...

>known_name_1

CACAGTCGAACGGC...
CACAGTCGAACGGC...
...

>unknown_name_2
...
```

"conv.csv" is a CSV-File with some more information:
```
unknown_name_1,unknown_dna_sequence_1,found_name_1,complete_dna_1,accession,own_score,bit-score,score,identity,hit-from,hit-to
unknown_name_2,unknown_dna_sequence_2,found_name_2,complete_dna_2,accession,own_score,bit-score,score,identity,hit-from,hit-to
...
```


## Settings
The default settings are listed in the "settings.json".
You can use your own settings by changing them or including a own json file like:

```console
user@pc:~$ python converter.py --input data/example.txt -o dna1 -c my_own_settings.json
```

## Help
Show help by typing:
```console
user@pc:~$ python converter.py --help
```

"settings.json":
```json
{
    "database":"nt",
    "default_score":"bit-score",
    "avoid_words":[
        "[Cc]omplete genome",
        "[Uu]ncultured",
        "[Uu]nknown"
    ],
    "score_words":{
        "16S.*RNA":1.5
    }
}
```

You can defined words or regex paterns in "avoid_words". Regex patterns or even words can be used to score occurence of them in the results by changing or extending "score_words". 

"default_score" can be set as:
<ul>
  <li>bit-score</li>
  <li>score</li>
  <li>identity</li>
  <li>hit-len</li>
</ul>

