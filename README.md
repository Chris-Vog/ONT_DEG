# Differential gene expression (**DEG**) analysis using Oxford Nanopore Technologies

Dies ist ein Workflow um differenziell exprimierte Gene mit Hilfe der Oxford Nanopore Technologie zu analysieren. Bei dieser Technologie handelt es sich um eine Methode, bei der lange Nukleotid-Sequenzen mithilfe von Nanoporen analysiert werden. Diese sind in einer Membran eingebettet und wo sie von einem Elektrolyt umgeben sind. Durch das Anlegen einer Spannung werden die Ionen in Bewegung gesetzt und folglich durch die Poren str�men. Die Bewegung der Ionen erzeugt einen stetigen Ionen-Fluss. Sobald DNA-Molek�le sich einer Pore n�hern, werden die doppelstr�ngigen Konstrukte aufgetrennt und aufgrund biophysischer Kr�fte durch die Pore geschleust. Der dadurch resultierende eingeschränkte Raum innerhalb der Pore f�hrt dazu, dass der Ionen-Fluss eingeschr�nkt ist und folglich reduziert wird. Diese Schwankungen hinsichtlich des Ionen-Flusses werden aufgenommen und k�nnen mithilfe von *Machine learning*-Algorithmen in Nukleotidsequenzen �bersetzt werden.

## Basecalling (Linux)

Unter dem Prozess des *Basecallings* versteht man die �bersetzung der Ver�nderungen im Stromfluss in eine Nukleotidsequenz. Die *Software* **Guppy** enth�lt den **ONT** Algorithmus sowie zus�tzliche Funktionen f�r die weitere prozessierung der Daten. Die Informationen hinsichtlich der Spannungs�nderung sind in sog. `.fast5`-Dateien gespeichert und die Sequenzen werden dann in `.fastq`-Dateien umgeschrieben. Weitere wichtige und n�tzliche Funktionen sind:

- Demultiplexing
- Adapter Trimming
- Qualit�tskontrolle

Da der Algorithmus sich zudem st�ndig weiterentwickelt, ist es ratsam m�glichst aktuelle Versionen zu verwenden. 

Die *Software* selbst l�sst sich mit drei essenziellen Befehlen ausf�hren.

|Funktion|Befehl|Input|
|----|-----|-----|
|Config-Datei f�r die Flow cell | -c / --config | dna_r9.4.1_450bps_sup.cfg|
|Pfad zu den .fast5-Dateien | -i / --input | ~path/to/files|
|Speicherort f�r die .fastq-Dateien | -s / --save_path | ~path/to/folder|

Informationen �ber s�mtliche Funktionen von Guppy sind in der Dokumentation aufgef�hrt.

    $ guppy_basecaller --help

Um den Prozess zu beschleunigen empfiehlt es sich diesen im GPU-Modus auszuf�hren. Hierf�r muss das Attribut `-x` definiert werden. Hierzu sollte der Befehl `auto` ausreichend sein. 

Der komplette Befehl:

    $ guppy_basecaller --min_qscore 7 --trim_barcodes --barcode_kits "SQK-PCB109" -i fast5/ -s fastq/ -c dna_r9.4.1_450bps_sup.cfg -x auto --gpu_runners_per_device 4 --chunks_per_runner 128

Als Resultat erh�lt man f�r jeden *Barcode* einen Ordner, der eine Vielzahl an `.fastq`-Dateien enth�lt. F�r die weitere Verarbeitung der Daten m�ssen diese konkateniert werden:

    $ cat path/to/fastq/{barcode}/*.fastq > filename.fastq

Dieser Befehl muss f�r jede Probe durchgef�hrt werden, sodass man f�r jede Probe eine einzige `.fastq`-Datei erh�lt, die alle Sequenzen enth�lt.

## Mapping using Minimap2 (Linux)

Beim `Mapping` werden kleinere Sequenzen mit einer langen Referenz-Sequenz abgeglichen und Bereichen zugeordnet. `Minimap2` ist dabei ein vielseitiges Programm, was insbesondere beim *mapping* langer Sequenzen von nicht allzu hoher Qualit�t gute Ergebnisse erzielt. 


