# Differential gene expression (**DEG**) analysis using Oxford Nanopore Technologies

Dies ist ein Workflow um differenziell exprimierte Gene mit Hilfe der Oxford Nanopore Technologie zu analysieren. Bei dieser Technologie handelt es sich um eine Methode, bei der lange Nukleotid-Sequenzen mithilfe von Nanoporen analysiert werden. Diese sind in einer Membran eingebettet und wo sie von einem Elektrolyt umgeben sind. Durch das Anlegen einer Spannung werden die Ionen in Bewegung gesetzt und folglich durch die Poren str?men. Die Bewegung der Ionen erzeugt einen stetigen Ionen-Fluss. Sobald DNA-Molek?le sich einer Pore n?hern, werden die doppelstr?ngigen Konstrukte aufgetrennt und aufgrund biophysischer Kr?fte durch die Pore geschleust. Der dadurch resultierende eingeschränkte Raum innerhalb der Pore f?hrt dazu, dass der Ionen-Fluss eingeschr?nkt ist und folglich reduziert wird. Diese Schwankungen hinsichtlich des Ionen-Flusses werden aufgenommen und k?nnen mithilfe von *Machine learning*-Algorithmen in Nukleotidsequenzen ?bersetzt werden.

## Bevor es losgeht

Die ben?tigten *Tools* wurden in einer Arbeitsumgebung installiert, die vor der Analyse aktiviert werden muss.

    $ source activate RNA-Seq

Des Weiteren wird f?r die Auswertung und insbesondere f?r die bioinformatische Aufbereitung der Daten eine bestimmte Ordnerstruktur benutzt. Diese kann zwar angepasst werden, aber die Nutzung dieser Ordnerstruktur vereinfacht den *Workflow*. Um die Ordnerstruktur zu erhalten muss folgender Befehl eurem Projektordner ausgef?hrt werden. Dazu macht ihr in eurem Projektordner einen Rechtsklick und w?hlt **Open in terminal** aus.

    $ mkdir -p ./{Analysis/{flagstat,Minimap2},ReferenceData,fast5,fastq}

F?r die weitere Analyse empfiehlt sich zudem die Befehle im Projektordner auszuf?hren!

## Basecalling (Linux)

Unter dem Prozess des *Basecallings* versteht man die ?bersetzung der Ver?nderungen im Stromfluss in eine Nukleotidsequenz. Die *Software* **Guppy** enth?lt den **ONT** Algorithmus sowie zus?tzliche Funktionen f?r die weitere prozessierung der Daten. Die Informationen hinsichtlich der Spannungs?nderung sind in sog. `.fast5`-Dateien gespeichert und die Sequenzen werden dann in `.fastq`-Dateien umgeschrieben. Weitere wichtige und n?tzliche Funktionen sind:

- Demultiplexing
- Adapter Trimming
- Qualit?tskontrolle

Da der Algorithmus sich zudem st?ndig weiterentwickelt, ist es ratsam m?glichst aktuelle Versionen zu verwenden. 

Die *Software* selbst l?sst sich mit drei essenziellen Befehlen ausf?hren.

|Funktion|Befehl|Input|
|----|-----|-----|
|Config-Datei f?r die Flow cell | -c / --config | dna_r9.4.1_450bps_sup.cfg|
|Pfad zu den .fast5-Dateien | -i / --input | ~path/to/files|
|Speicherort f?r die .fastq-Dateien | -s / --save_path | ~path/to/folder|

Informationen ?ber s?mtliche Funktionen von Guppy sind in der Dokumentation aufgef?hrt.

    $ guppy_basecaller --help

Um den Prozess zu beschleunigen empfiehlt es sich diesen im GPU-Modus auszuf?hren. Hierf?r muss das Attribut `-x` mit dem Befehl `auto` ausgef?hrt werden. 

Der komplette Befehl:

    $ guppy_basecaller --min_qscore 7 --trim_barcodes --barcode_kits "SQK-PCB109" -i fast5/ -s fastq/ -c dna_r9.4.1_450bps_sup.cfg -x auto --gpu_runners_per_device 4 --chunks_per_runner 128

Als Resultat erh?lt man f?r jeden *Barcode* einen Ordner, der eine Vielzahl an `.fastq`-Dateien enth?lt. F?r die weitere Verarbeitung der Daten m?ssen diese konkateniert werden:

    $ cat path/to/fastq/{barcode}/*.fastq > filename.fastq

Dieser Befehl muss f?r jede Probe durchgef?hrt werden, sodass man f?r jede Probe eine einzige `.fastq`-Datei, die alle Sequenzen enth?lt.

Um zudem Platz zu sparen k?nnen die Dateien komprimiert werden. Die hier verwendeten *Tools* sind aber in der Lage sowohl komprimierte als auch nicht komprimierte Dateien zu analysieren.

    $ gzip path/to/.fastq-File

### Quality-Control mittels MinIONQC (R auf verschiedenen OS anwendbar)

`MinIONQC` ist ein ausgearbeitetes **R**-Skript zur Analyse von ONT-Datens?tzen. Hierzu wird lediglich die `sequencing_summary.txt`-Datei ben?tigt. Mit dem folgenden Befehl wird `MinIONQC` ausgef?hrt:

    $ Rscript MinIONQC.R -i path/to/parent/directory -q 7

Erl?uterungen:
- **MinIONQC**: Pfad zum Skript
- **path/to/parent/directory**: Pfad zum Input-Ordner, der die `sequencing_summary.txt`-Datei enth?lt.
- **-q/ --qscore_cutoff**: Grenzwert f?r den durchschnittlichen *Mean Q Score* (7 ist hier voreingestellt)


## Mapping using Minimap2 (Linux)

Beim `Mapping` werden kleinere Sequenzen mit einer langen Referenz-Sequenz abgeglichen und Bereichen zugeordnet. `Minimap2` ist dabei ein vielseitiges Programm, was insbesondere beim *mapping* langer Sequenzen von nicht allzu hoher Qualit?t gute Ergebnisse erzielt. 

### Download des Referenz-Genoms und der Annotation (Linux)

Die Sequenzier-Daten werden bei diesem Schritt mit dem humanen Genom abgeglichen und die Position mit der h?chsten ?bereinstimmung identifiziert. Das humane Genom wird als `FASTA` zur Verf?gung gestellt.

    $ wget -P path/to/ReferenceData https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.p13.genome.fa.gz

Nukleotid-Sequenz der GRCh38.p13 Genomassemblierung aller Regionen, inklusive Referenz-Chromosomen, Korrekturen und Haplotypen. Die dazugeh?rige Annotation ist in einer `.gtf`-Datei gespeichert. Wichtig hierbei ist, dass die Versionen ?bereinstimmen.

    $ wget -P path/to/ReferenceData https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.chr_patch_hapl_scaff.annotation.gtf.gz

### Indexing (Linux)

Ohne Angaben von weiteren Informationen verwendet `Minimap2` eine Referenz-Datenbank und eine Sequenz-Datei und f?hrt ein ann?herndes *Mapping* durch. F?r das humane Genom muss jedoch vorerst ein Index generiert werden.

Verwendetes Referenz-Genom:

    Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

Es empfiehlt sich immer eine definierte Ordnerstruktur einzuhalten. S?mtliche Referenzdaten, befinden sich hier im Ordner `ReferenceData`. Der Index wird folgenderma?en erstellt.:

    $ minimap2 -t 16 -k14 -w5 -d ReferenceData/reference.mmi ReferenceData/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

### Map long noisy reads (Linux)

F?r alle Anwendungen verwendet `Minimap2` den selben Algorithmus. Dennoch muss der Befehl in Abh?ngigkeit der verwendeten Sequencing-Technologie angepasst werden:

### Post-processing of sequencing data using Samtools (Linux)

`Samtools` ist eine Sammplung von Programmen, um mit Sequenzierdaten zu interagieren und weiterzubearbeiten. Es wird hier genutzt die Sequenzen in Abh?ngigkeit ihrer Lokalisation im Genom zu sortieren:

    $ minimap2 -t {threads} -ax splice -k14 --secondary=no path/to/reference.mmi {input.fastq} | samtools view -Sb | samtools sort -o {output.sorted.bam}

Mithilfe des *Pipe Operators* `|` k?nnen die Ergebnisse einer Operation weiteren Befehlen ?bergeben werden. 

Der {input} und der {output} sollten dabei identisch sein.

Die Qualit?t und die Effizienz des Mappings wird ebenfalls mittels samtools analysiert. Hierzu werden die flagstats erstellt, die die Anzahl erfolgreicher Alignments aufz?hlen.

    $ samtools flagstat /path/to/bam_file > {flagstat}.txt

Die so erstellten .txt-Dateien werden mit Hilfe von **R** ausgewertet und visualisiert.

Des Weiteren k?nnen *Mappings* mit `Samtools` indexiert werden. Dies ist notwendig, wenn die Daten beispielsweise mit dem **IGV browser** analysiert werden. 

    $ samtools index -b /path/to/{input}.sorted.bam /path/to/{input}.sorted.bam.bai

Die resultierende Datei ist noch nicht vorhanden und wir mit Hilde diesen Befehls benannt und erstellt. Verwende hierzu den gleichen Dateinamen und f?ge `.bai` hinzu. 