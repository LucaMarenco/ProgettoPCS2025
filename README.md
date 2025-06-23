Progetto di Programmazione e Calcolo Scientifico
Luca Marenco, Alessandro Contini, Massimo Cirronis
2025

# 1 Poliedri Geodetici e i loro duali
Il progetto `e costituito da due parti. Dato in input una quadrupla di numeri
interi (p, q, b, c):
1. Definire i solidi geodetici di classe I e i loro duali.
2. Definire i solidi geodetici di classe II e i loro duali.

## 1.1 Parte I
La prima parte del progetto consiste nel:

1. Definire una struttura dati che permetta la memorizzazione di tutte le
propriet`a di un poliedro. Ogni poliedro deve essere rappresentato dalle:

• Celle 0D oppure vertici. Ogni vertice `e caratterizzato da un identificativo
(numero sequenziale a partire da 0), le 3 coordinate (x,y,z ).

• Cella 1D oppure lati. Ciascun lato `e caratterizzato da un identificativo
(numero sequenziale a partire da 0) e dagli ID dei vertici
di origine e fine che individuano in maniera univoca il lato nel
poliedro.

• Cella 2D oppure faccia. Ciascuna faccia `e caratterizzata da un
identificativo (numero sequenziale a partire da 0), dal numero di
vertici e di lati, e da due liste definite dagli ID dei vertici e dei
lati che individuano in maniera univoca la faccia nel poliedro. In
ciascuna faccia, i lati e i vertici devono essere ordinati in modo tale
che, a meno dell’orientamento del lato, risulti
faces.edges[e].end == faces.edges[(e + 1)%E].origin
faces.vertices[e] == faces.edges[e].origin

• Cella 3D oppure poliedro. Ciascun poliedro `e caratterizzato da un
identificativo (numero sequenziale a partire da 0), dal numero di
vertici, di lati e di facce, e da tre liste definite dagli ID dei vertici,
dei lati e delle facce che individuano in maniera univoca il poliedro.

2. Dato in input una quadrupla di numeri interi (p,q,b,0) oppure (p,q,0,b),
se p = 3 il programma deve restituire in output il poliedro geodetico
di classe I corrispondente. In particolare, il programma deve restituire
4 file .txt denominati Cell0Ds.txt, Cell1Ds.txt, Cell2Ds.txt e
Cell3Ds.txt riportanti le principali propriet`a che caratterizzano le varie
celle del poliedro. Inoltre, il programma deve consentire la stampa su
Paraview dei vertici e dei lati del poliedro. Non `e richiesta la stampa
delle facce e del poliedro stesso.

3. Dato in input una quadrupla di numeri interi (p,q,b,0) oppure (p,q,0,b),
se q = 3 il programma deve restituire in output il poliedro di Goldberg
di classe I corrispondente, consentendo gli output previsti dal punto
precedente.

4. Costruire tali poliedri in modo tale che tutti i loro vertici giacciano sulla
sfera di raggio 1 centrata nell’origine degli assi cartesiani.

5. Inoltre, se l’input `e definito da una 6-tupla di numeri
(p,q,b,c,id vertice 1,id vertice 2),
trovare un cammino minimo che unisce i vertici contrassegnati dagli identificativi
id vertice 1 e id vertice 2 (se validi) sul grafo avente come nodi
i vertici del poliedro e come lati le celle 1D dello stesso. Il cammino minimo
andr`a evidenziato su Paraview, assegnando propriet`a ShortPath = 1
ai lati e ai vertici che appartengono al cammino minimo e ShortPath = 0
a quelli che non vi appartengono. Inoltre, il programma deve stampare a
schermo il numero di lati che compongono il cammino minimo e la somma
delle loro lunghezze.

6. Verificare sempre la correttezza dell’input.

Per ogni unit`a logica `e necessario verificarne il corretto funzionamento utilizzando
i GoogleTest

## 1.2 Parte II
La seconda parte del progetto consiste nel modificare l’algoritmo prodotto nella
parte I del progetto, in modo tale da consentire la costruzione dei solidi geodetici
di classe II (b = c).
