# FASE 15 — Confronto diagnostico Cp a quadratura comune

**Solo postprocessing, nessuna nuova simulazione CFD.** Il campo primitivo è
quello native/PROJECTED, CFL=.05, dopo 9000 iterazioni della FASE 15D:
**non è certificato pienamente stazionario**. Il riferimento è il campo
convergente EULER nativo della FASE 15A, non una soluzione analitica esatta.

## Risposte alle tre domande

1. **Cp coincide ragionevolmente lontano dal TE?** Buona somiglianza nella
   parte centrale .10≤x/c≤.90: RMS campionaria di DeltaCp circa .0047 e
   massimo circa .014 su ciascun ramo. Non è però un accordo uniforme su
   tutto il tratto lontano dal TE: al LE la RMS è circa .035 e il massimo
   .070. Sull'intero x/c≤.90 la RMS campionaria è .01844 (pesata sulla
   lunghezza .01068). La sovrapposizione visiva globale non esclude uno
   scostamento sistematico di pressione al naso.
2. **La differenza di CD è concentrata al TE?** No: **84.34% viene dal LE**,
   12.42% dal TE e 3.25% dalla parte centrale. Il massimo locale di DeltaCp
   è al TE, ma non domina l'integrale assiale. Le due mezze facce possedute
   dal nodo TE danno addirittura un piccolo contributo positivo al DeltaCD,
   opposto al suo valore totale negativo.
3. **Il confronto aerodinamico è soddisfacente senza dichiarare convergenza?**
   Solo parzialmente: la forma della distribuzione e il Cp centrale sono
   simili, ma la resistenza non è soddisfacente e le differenze LE/TE restano
   significative. Non si valida l'aerodinamica globale né lo stato stazionario.
   Somiglianza dei campi, accuratezza di CD e convergenza del residuo sono
   proprietà distinte.

## Dati e significato delle pressioni

Percorsi effettivamente verificati, relativi al repository:

- `QuickStart/mesh_NACA0012_inv.su2`;
- `TestCases/primitive_euler_examples/naca0012_aoa0/results/native/native_state_000599.dat`;
- `.../results15d/native_initial.dat`, conversione Float64 della FASE 15D;
- `.../results15d/run_final_native_PROJECTED_cfl005/primitive_state_00000.dat`
  e `primitive_state_09000.dat`;
- nella stessa cartella, `surface_primitive_pressure_0.csv`,
  `surface_primitive_pressure_9000.csv` e `history.csv`.

I due CSV di superficie sono scritti con 17 cifre significative da
`CPrimitiveEulerOutput::WriteAdditionalFiles`, che chiama
`CPrimitiveEulerSolver::ComputeAerodynamics`. Il primo è già l'estrazione
del campo nativo importato **prima del primo aggiornamento**; il secondo
è l'estrazione dello stato a 9000. Entrambi usano la stessa ricostruzione
primitiva WLS/MUSCL affine, limiter NONE, gli stessi punti delle mezze facce
e le stesse normali integrate. Non occorre rieseguire il solver o ricostruire
le pressioni con una nuova implementazione Python.

La pressione del restart nativo è stata verificata con
`p=(gamma-1)*(rhoE-|rho*u|²/(2*rho))`, gamma=1.4. Density e Pressure
iniziali importate coincidono con la conversione Float64. Le piccole
proiezioni iniziali delle velocità ordinarie (~1e-10, documentate in 15D)
non cambiano rho/p. La dicitura **nativo a quadratura comune** indica questo
campo importato, valutato tramite la ricostruzione primitiva; non indica
il Cp nodale prodotto dal postprocessing nativo.

La fonte principale non è un VTU Float32. Si usano i CSV di pressione
Float64 e i restart Float64; i campi archiviati sono letti e mai modificati.
Normalizzazione: p_inf=1, q_inf=.028, REF_AREA=1, AOA=0°.

## Geometria, accoppiamento e controlli

400 mezze facce del marker `airfoil`, 200 per dorso/ventre; 200 vertici di
parete, ciascuno proprietario di due campioni. Per ogni lato primale a–b,
i punti indipendenti attesi sono (3*a+b)/4 e (a+3*b)/4. Le normali sono
ricostruite geometricamente dal lato, con modulo metà della lunghezza e
verso opposto al triangolo fluido adiacente. Tutti i lati hanno esattamente
un triangolo fluido incidente.

Le normali così ricavate coincidono con quelle nei CSV; le coordinate
ricalcolate differiscono al massimo per 1.11e-16 in x e 6.94e-18 in y.
**I due CSV fra loro hanno coordinate, ID e normali identici**, consentendo
un accoppiamento uno a uno senza interpolazione. Nel CSV finale sono
conservate le effettive coordinate comuni archiviate. La somma delle
normali chiude entro roundoff. Tutti i valori sono finiti e le pressioni
positive; le formule Cp/pressione sono verificate.

Dorso: y_b>0; ventre: y_b<0. Nessun punto di quadratura ha y_b=0. I due rami
sono ordinati separatamente per x/c, monotono su questa mesh NACA. Non si
collegano i rami fra loro. Il nodo TE 199 possiede due campioni distinti,
uno per ramo: nessun doppio conteggio. Linee e marker collegano soltanto
campioni effettivi; nessuno smoothing o interpolazione nei grafici o nelle
metriche/integrali.

## Metriche Cp

DeltaCp=Cp_primitivo−Cp_nativo. Media e RMS campionarie assegnano lo stesso
peso a ogni mezza faccia; le versioni pesate usano
ell_b=sqrt(S_x²+S_y²), cioè la lunghezza reale della mezza faccia:

    mean_sample = sum DeltaCp / n
    RMS_sample = sqrt(sum DeltaCp² / n)
    mean_length = sum ell_b*DeltaCp / sum ell_b
    RMS_length = sqrt(sum ell_b*DeltaCp² / sum ell_b)

| Ramo | Regione | n | Media campionaria | RMS campionaria | Max assoluto | x/c del max | RMS su lunghezza |
|---|---|---:|---:|---:|---:|---:|---:|
| Dorso | LE, x<.10 | 42 | -.0305696 | .0352424 | .0697621 | .0000150 | .0278249 |
| Dorso | .10≤x≤.90 | 117 | -.00342659 | .00468694 | .0140350 | .103560 | .00437784 |
| Dorso | TE, x>.90 | 41 | .0465388 | .0957953 | .319804 | .998690 | .0484970 |
| Ventre | LE, x<.10 | 42 | -.0299392 | .0347808 | .0697128 | .0000150 | .0273021 |
| Ventre | .10≤x≤.90 | 117 | -.00347356 | .00475296 | .0139839 | .103560 | .00445168 |
| Ventre | TE, x>.90 | 41 | .0466437 | .0958565 | .320002 | .998690 | .0485534 |

Su tutto il dorso: media .00111629, RMS .0464210, RMS su lunghezza .0184264.
Sul ventre: media .00124268, RMS .0463777, RMS su lunghezza .0183652.
Tutte le medie pesate, le posizioni y dei massimi e i dati aggregati sono
in [cp_metrics.csv](cp_metrics.csv). La differenza fra RMS campionaria e
pesata risente dell'addensamento della mesh al TE: non sono norme intercambiabili.

La differenza centrale è piccola ma distribuita; al LE è prevalentemente
negativa e sistematica, non un singolo outlier. Vicino al TE compaiono scarti
più grandi e variazioni locali, resi visibili dagli zoom con marker.
Il riferimento ricostruito non è una previsione analitica di Cp: eventuali
picchi comuni ai due campi non sono automaticamente validati fisicamente.

## Forze e origine del CD negativo

Integrazione diretta senza interpolazione grafica:

    Fx = sum_b (p_b-1)*S_x,b
    CD = Fx/.028
    DeltaCD_G = sum_{b in G} (p_primitive,b-p_native,b)*S_x,b/.028

| Campo, quadratura comune | Fx | CD |
|---|---:|---:|
| Nativo importato, pressione ricostruita | +8.545041877e-6 | +.000305180067034459 |
| Primitivo a 9000, pressione ricostruita | -8.998151167e-5 | -.003213625416941096 |
| Differenza | -9.852655355e-5 | **-.003518805483975554** |

Massimo scarto del CD ricalcolato dalla history primitiva iniziale/finale:
**5.10e-18**. Anche CL è riprodotto entro roundoff. Le regioni ricostruiscono
il totale entro la precisione numerica; percentuali firmate rispetto al
DeltaCD totale, qui ben diverso da zero:

| Regione | n | Fx nativo | Fx primitivo | DeltaFx | DeltaCD | Quota firmata |
|---|---:|---:|---:|---:|---:|---:|
| LE, x<.10 | 84 | +1.110987859e-4 | +2.800322717e-5 | -8.309555877e-5 | -.002967698527 | 84.3382% |
| Centrale, .10≤x≤.90 | 234 | +1.554529658e-5 | +1.234685630e-5 | -3.198440273e-6 | -.000114230010 | 3.2463% |
| TE, x>.90 | 82 | -1.180990406e-4 | -1.303315951e-4 | -1.223255451e-5 | -.000436876947 | 12.4155% |

Il CD nativo comune è il piccolo resto di contributi anteriori e posteriori
quasi compensati. La perdita della forza assiale positiva anteriore domina
la variazione: uno scarto Cp di poche centesime al LE può quindi compromettere
CD anche quando le curve globali appaiono vicine. Norme di DeltaCp, massimo
locale e integrale con S_x misurano proprietà differenti.

### Ultimi due segmenti del TE

Sono mostrate tutte e quattro le mezze facce dei due segmenti adiacenti al
TE, includendo i due proprietari vicini e i due campioni del nodo 199:

| Ramo | Proprietario | x/c | DeltaCp | Contributo a DeltaCD |
|---|---:|---:|---:|---:|
| Ventre | 0 | .999812514 | +.1001201 | -1.818629812e-6 |
| Ventre | TE 199 | .999937505 | -.02236237 | +4.062008611e-7 |
| Dorso | 198 | .999812514 | +.1003723 | -1.823210925e-6 |
| Dorso | TE 199 | .999937505 | -.02236322 | +4.062163629e-7 |

Somma dei quattro contributi: circa -2.82942e-6, solo .0804% del DeltaCD
totale. I due del solo nodo TE sommano +8.12417e-7: contrastano leggermente
il DeltaCD negativo. Il maggior contributo negativo della regione TE è
in .99<x≤.999 (-.000315658), non nelle due ultime mezze facce del nodo TE.
Dettagli in [te_adjacent_half_faces.csv](te_adjacent_half_faces.csv) e
[trailing_edge_bins.csv](trailing_edge_bins.csv).

### Verifica indipendente della convenzione nativa

`CFVMFlowSolverBase::Pressure_Forces`, in
`SU2_CFD/include/solvers/CFVMFlowSolverBase.inl:1922–1958`, usa pressione
**nodale** e normali dei vertex, non pressione ricostruita ai punti delle
mezze facce. Applicare p_nodale del restart nativo alle medesime mezze facce
e sommarle per proprietario riproduce la convenzione dei vertex, con il
segno coerente delle normali.

Questo controllo Float64 riproduce CD_nativo=.000409850480820 circa e il CL
archiviato, entro 3e-14. La differenza dal suo CD a ricostruzione comune
(.000305180067034) deriva dal diverso campionamento della pressione; non
è stata forzata una coincidenza fra definizioni diverse. Il confronto
principale usa sempre pressione **ricostruita comune** per entrambi i campi.
La piccolezza del CD nativo non rende il solver un riferimento infallibile
di accuratezza fisica: qui è un controllo di coerenza e un confronto di campi.

## Limiti e artefatti

Il residuo finale resta quello della 15D: globale 2.77177e-5, TE .0794207,
con variazioni residue di CL/CD. L'origine spaziale del DeltaCD è qui
localizzata, ma non si dimostra se la differenza sia causata da rilassamento
incompleto, formulazione primitiva, discretizzazione o politica TE. Una
concentrazione del **residuo** vicino al TE non implicava una concentrazione
della **differenza di forza** nello stesso luogo: i due dati non si contraddicono.

Il precedente `results15d/cp_comparison_status.txt` resta False e invariato.
Questa estrazione è stata esplicitamente richiesta come diagnostica, non
come superamento del gate di stazionarietà.

- [Grafico principale PNG](cp_native_vs_primitive.png) e
  [PDF vettoriale](cp_native_vs_primitive.pdf), 4 curve su rami separati;
- [zoom TE .90–1](cp_trailing_edge_zoom.png),
  [dettaglio TE .98–1](cp_trailing_edge_detail.png) e
  [zoom LE 0–.10](cp_leading_edge_zoom.png), anche PDF;
- [DeltaCp](delta_cp.png), [contributi regionali CD](drag_difference_by_region.png)
  e [contributo cumulativo](drag_difference_cumulative.png);
- `cp_comparison.csv`: 400 coppie comuni, pressioni, Cp, normali e contributi
  di forza; `drag_decomposition.csv`, `forces.csv`, metriche e dettagli TE;
- `provenance.json`: sorgenti e SHA256, precisione, normalizzazione e
  controlli. Gli hash degli input sono ricontrollati a fine estrazione.

Riproduzione, senza processi CFD:

```sh
python3 UnitTests/postprocess_primitive_euler_cp15.py --output-dir /private/tmp/nuovo-confronto-cp
python3 -m py_compile UnitTests/postprocess_primitive_euler_cp15.py
git diff --check
```

La directory di uscita deve essere nuova: il comando rifiuta di sovrascrivere
artefatti esistenti. Eseguita anche una ripetizione indipendente in directory
temporanea; controlli di geometria, campionamento, conversione e forze
passati. Nessuna modifica produttiva, nessuna simulazione, nessuna
alterazione dei campi o degli artefatti 15A–15D precedenti, nessuna operazione
Git di scrittura. Lavoro precedente integralmente conservato nel working tree.
