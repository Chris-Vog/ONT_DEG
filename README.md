# Differential gene expression (**DEG**) analysis using Oxford Nanopore Technologies

Dies ist ein Workflow um differenziell exprimierte Gene mit Hilfe der Oxford Nanopore Technologie zu analysieren. Bei dieser Technologie handelt es sich um eine Methode, bei der lange Nukleotid-Sequenzen mithilfe von Nanoporen analysiert werden. Diese sind in einer Membran eingebettet und wo sie von einem Elektrolyt umgeben sind. Durch das Anlegen einer Spannung werden die Ionen in Bewegung gesetzt und folglich durch die Poren strÃ¶men. Die Bewegung der Ionen erzeugt einen stetigen Ionen-Fluss. Sobald DNA-MolekÃ¼le sich einer Pore nÃ¤hern, werden die doppelstrÃ¤ngigen Konstrukte aufgetrennt und aufgrund biophysischer KrÃ¤fte durch die Pore geschleust. Der dadurch resultierende eingeschrÃ¤nkte Raum innerhalb der Pore fÃ¼hrt dazu, dass der Ionen-Fluss eingeschrÃ¤nkt ist und folglich reduziert wird. Diese Schwankungen hinsichtlich des Ionen-Flusses werden aufgenommen und kÃ¶nnen mithilfe von *Machine learning*-Algorithmen in Nukleotidsequenzen Ã¼bersetzt werden.

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

