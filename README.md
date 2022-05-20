# Differential gene expression (**DEG**) analysis using Oxford Nanopore Technologies

Dies ist ein Workflow um differenziell exprimierte Gene mit Hilfe der Oxford Nanopore Technologie zu analysieren. Bei dieser Technologie handelt es sich um eine Methode, bei der lange Nukleotid-Sequenzen mithilfe von Nanoporen analysiert werden. Diese sind in einer Membran eingebettet und wo sie von einem Elektrolyt umgeben sind. Durch das Anlegen einer Spannung werden die Ionen in Bewegung gesetzt und folglich durch die Poren strömen. Die Bewegung der Ionen erzeugt einen stetigen Ionen-Fluss. Sobald DNA-Moleküle sich einer Pore nähern, werden die doppelsträngigen Konstrukte aufgetrennt und aufgrund biophysischer Kräfte durch die Pore geschleust. Der dadurch resultierende eingeschrÃ¤nkte Raum innerhalb der Pore führt dazu, dass der Ionen-Fluss eingeschränkt ist und folglich reduziert wird. Diese Schwankungen hinsichtlich des Ionen-Flusses werden aufgenommen und können mithilfe von *Machine learning*-Algorithmen in Nukleotidsequenzen übersetzt werden.

## Basecalling (Linux)

Unter dem Prozess des *Basecallings* versteht man die Übersetzung der Veränderungen im Stromfluss in eine Nukleotidsequenz. Die *Software* **Guppy** enthält den **ONT** Algorithmus sowie zusätzliche Funktionen für die weitere prozessierung der Daten. Die Informationen hinsichtlich der Spannungsänderung sind in sog. `.fast5`-Dateien gespeichert und die Sequenzen werden dann in `.fastq`-Dateien umgeschrieben. Weitere wichtige und nützliche Funktionen sind:

- Demultiplexing
- Adapter Trimming
- Qualitätskontrolle

Da der Algorithmus sich zudem ständig weiterentwickelt, ist es ratsam möglichst aktuelle Versionen zu verwenden. 

Die *Software* selbst lässt sich mit drei essenziellen Befehlen ausführen.

|Funktion|Befehl|Input|
|----|-----|-----|
|Config-Datei für die Flow cell | -c / --config | dna_r9.4.1_450bps_sup.cfg|
|Pfad zu den .fast5-Dateien | -i / --input | ~path/to/files|
|Speicherort für die .fastq-Dateien | -s / --save_path | ~path/to/folder|

Informationen über sämtliche Funktionen von Guppy sind in der Dokumentation aufgeführt.

    $ guppy_basecaller --help

Um den Prozess zu beschleunigen empfiehlt es sich diesen im GPU-Modus auszuführen. Hierfür muss das Attribut `-x` mit dem Befehl `auto` ausgeführt werden. 

Der komplette Befehl:

    $ guppy_basecaller --min_qscore 7 --trim_barcodes --barcode_kits "SQK-PCB109" -i fast5/ -s fastq/ -c dna_r9.4.1_450bps_sup.cfg -x auto --gpu_runners_per_device 4 --chunks_per_runner 128

Als Resultat erhält man für jeden *Barcode* einen Ordner, der eine Vielzahl an `.fastq`-Dateien enthält. Für die weitere Verarbeitung der Daten müssen diese konkateniert werden:

    $ cat path/to/fastq/{barcode}/*.fastq > filename.fastq

Dieser Befehl muss für jede Probe durchgeführt werden, sodass man für jede Probe eine einzige `.fastq`-Datei erhält, die alle Sequenzen enthält.

### Quality-Control mittels MinIONQC (R auf verschiedenen OS anwendbar)

`MinIONQC` ist ein ausgearbeitetes **R**-Skript zur Analyse von ONT-Datensätzen. Hierzu wird lediglich die `sequencing_summary.txt`-Datei benötigt. Mit dem folgenden Befehl wird `MinIONQC` ausgeführt:

    $ Rscript MinIONQC.R -i path/to/parent/directory -q 7

Erläuterungen:
- **MinIONQC**: Pfad zum Skript
- **path/to/parent/directory**: Pfad zum Input-Ordner, der die `sequencing_summary.txt`-Datei enthält.
- **-q/ --qscore_cutoff**: Grenzwert für den durchschnittlichen *Mean Q Score* (7 ist hier voreingestellt)


## Mapping using Minimap2 (Linux)

Beim `Mapping` werden kleinere Sequenzen mit einer langen Referenz-Sequenz abgeglichen und Bereichen zugeordnet. `Minimap2` ist dabei ein vielseitiges Programm, was insbesondere beim *mapping* langer Sequenzen von nicht allzu hoher Qualität gute Ergebnisse erzielt. 

### Download des Referenz-Genoms und der Annotation (Linux)

Die Sequenzier-Daten werden bei diesem Schritt mit dem humanen Genom abgeglichen und die Position mit der höchsten Übereinstimmung identifiziert. Das humane Genom wird als `FASTA` zur Verfügung gestellt.

    $ wget -P path/to/ReferenceData https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.p13.genome.fa.gz

Nukleotid-Sequenz der GRCh38.p13 Genomassemblierung aller Regionen, inklusive Referenz-Chromosomen, Korrekturen und Haplotypen. Die dazugehörige Annotation ist in einer `.gtf`-Datei gespeichert. Wichtig hierbei ist, dass die Versionen übereinstimmen.

    $ wget -P path/to/ReferenceData https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.chr_patch_hapl_scaff.annotation.gtf.gz

### Indexing (Linux)

Ohne Angaben von weiteren Informationen verwendet `Minimap2` eine Referenz-Datenbank und eine Sequenz-Datei und führt ein annäherndes *Mapping* durch. Für das humane Genom muss jedoch vorerst ein Index generiert werden.

Verwendetes Referenz-Genom:

    Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

Es empfiehlt sich immer eine definierte Ordnerstruktur einzuhalten. Sämtliche Referenzdaten, befinden sich hier im Ordner `ReferenceData`. Der Index wird folgendermaßen erstellt.:

    $ minimap2 -t 16 -k14 -w5 -d ~/ReferenceData/reference.mmi ~/ReferenceData/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

### Map long noisy reads (Linux)

Für alle Anwendungen verwendet `Minimap2` den selben Algorithmus. Dennoch muss der Befehl in Abhängigkeit der verwendeten Sequencing-Technologie angepasst werden:

    $ minimap2 -t {threads} -ax splice -k14 --secondary=no {input.index} {input.fastq}

### Post-processing of sequencing data using Samtools (Linux)

`Samtools` ist eine Sammplung von Programmen, um mit Sequenzierdaten zu interagieren und weiterzubearbeiten. Es wird hier genutzt die Sequenzen in Abhängigkeit ihrer Lokalisation im Genom zu sortieren:

    $ minimap2 -t {threads} -ax splice -k14 –secondary=no {input.index} {input.fastq} | samtools view -Sb | samtools sort -o {output.sorted.bam}

Hier ist zu beachten, dass der vorherige Befehl sich hier wiederholt. Mithilfe des *Pipe Operators* `|` können die Ergebnisse einer Operation weiteren Befehlen übergeben werden.

Die Qualität und die Effizienz des Mappings wird ebenfalls mittels samtools analysiert. Hierzu werden die flagstats erstellt, die die Anzahl erfolgreicher Alignments aufzählen.

    $ samtools flagstat /path/to/bam_file > {flagstat}.txt

Die so erstellten .txt-Dateien werden mit Hilfe von **R** ausgewertet und visualisiert.

Des Weiteren können *Mappings* mit `Samtools` indexiert werden. Dies ist notwendig, wenn die Daten beispielsweise mit dem **IGV browser** analysiert werden. 

    $ samtools index -b /path/to/{input}.sorted.bam /path/to/{input}.sorted.bam.bai

Die resultierende Datei ist noch nicht vorhanden und wir mit Hilde diesen Befehls benannt und erstellt. Verwende hierzu den gleichen Dateinamen und füge `.bai` hinzu. 