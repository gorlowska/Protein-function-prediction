# Meeting 2 — raport postępu (Team 2)

## Cel
Przygotowanie danych do przewidywania funkcji białek z sekwencji: GO jako wiele etykiet.
Organizmy: Candida albicans, Candida parapsilosis, Saccharomyces cerevisiae.

## Źródła i odtwarzalność (około 45 s)
GAF UniProt-centric z GO i reviewed Swiss-Prot z UniProt. Zapisano raw,
daty UTC, wersje, zapytania i SHA256. Eksport JSON GO uzupełnia TSV o dowody.
Wydanie JSON musi odpowiadać TSV; zestawy accession i GO są sprawdzane.
Seed podziału: 42. Wersja i polecenia MMseqs2 są w osobnych raportach.

## Porównanie i filtr dowodów (około 90 s)
Porównano bezpośrednie pozytywne GO per białko w obu źródłach.
Pełne zestawy UniProt pozostają w master jako ground truth przed filtrem.
Końcowe etykiety są ich podzbiorem, z dozwolonym dowodem dla konkretnej pary
białko–GO w UniProt lub pasującym GAF. Dokładne konflikty GAF NOT wykluczono.
Nie przenoszono adnotacji innej izoformy. Wykluczono białka bez końcowej etykiety.
Nie usuwano białek wyłącznie dlatego, że nie mają odpowiednika w GAF.

```text
                 species  reviewed  with_matching_positive_GAF  identical_GO_sets  retained  excluded
        Candida albicans      1466                        1192                 83       704       762
    Candida parapsilosis        23                           7                  0         9        14
Saccharomyces cerevisiae      7919                        6057                845      5543      2376
```

## Homologia i podział (około 90 s)
FASTA wszystkich zachowanych sekwencji przetworzono przez MMseqs2 easy-cluster:
min-seq-id 0.5, coverage 0.8, cov-mode 0; alignment-mode 3, cluster-reassign 1.
Całe klastry scalono dodatkowo po identycznych sekwencjach i wykrytych parach
spełniających progi w all-vs-all easy-search. Nie dzielono tych grup.
Podział: 128 prób zachłannych z seed 42, cel 8:1:1 liczby białek,
uwzględnienie rozkładu gatunków i częstych etykiet.

Klastry MMseqs2: 5504; końcowe grupy: 5481.

```text
split  proteins  clusters  label_pairs  fraction
train      5005      4382        31647  0.800032
 eval       625       550         3791  0.099904
 test       626       549         3773  0.100064
```

## Walidacja, ograniczenia i kolejny meeting (około 75 s)
Kontrole zakończyły się poprawnie: brak rozdzielonych klastrów, brak identycznych
sekwencji między zbiorami, brak kwalifikujących trafień wyszukiwania między zbiorami,
kompletne FASTA i zgodność etykiet z polityką dowodów.
Wyszukiwanie MMseqs2 jest heurystyczne; nie jest wyczerpującym dowodem braku każdej podobnej pary.
Silna nierównowaga gatunków i mało danych C. parapsilosis ograniczają wnioskowanie.
Rzadkie terminy GO nie muszą występować we wszystkich zbiorach.
Słownik klas wyznaczono wyłącznie na train; pełne etykiety eval/test zachowano.
Brak adnotacji nie oznacza biologicznego braku funkcji.

GO nieobecne w train: 220 w eval i 236 w test.

Następnie meeting 3: działający szkielet wejście/wyjście na gotowych plikach.

Pomiar maksymalnego RSS procesów MMseqs2 (KiB, nie całego WSL):

```text
clusters: Maximum resident set size (kbytes): 905936
homology_audit: Maximum resident set size (kbytes): 928868
```
