# Differential gene expression (**DEG**) analysis using Oxford Nanopore Technologies

Dies ist ein Workflow um differenziell exprimierte Gene mit Hilfe der Oxford Nanopore Technologie zu analysieren. Bei dieser Technologie handelt es sich um eine Methode, bei der lange Nukleotid-Sequenzen mithilfe von Nanoporen analysiert werden. Diese sind in einer Membran eingebettet und wo sie von einem Elektrolyt umgeben sind. Durch das Anlegen einer Spannung werden die Ionen in Bewegung gesetzt und folglich durch die Poren strömen. Die Bewegung der Ionen erzeugt einen stetigen Ionen-Fluss. Sobald DNA-Moleküle sich einer Pore nähern, werden die doppelsträngigen Konstrukte aufgetrennt und aufgrund biophysischer Kräfte durch die Pore geschleust. Der dadurch resultierende eingeschränkte Raum innerhalb der Pore führt dazu, dass der Ionen-Fluss eingeschränkt ist und folglich reduziert wird. Diese Schwankungen hinsichtlich des Ionen-Flusses werden aufgenommen und können mithilfe von *Machine learning*-Algorithmen in Nukleotidsequenzen übersetzt werden.

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

