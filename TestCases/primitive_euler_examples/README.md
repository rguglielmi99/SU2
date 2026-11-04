# Primitive Euler — esempi SU2_CFD nativi

Questi casi usano il driver standard e mesh `.su2` già incluse. Non serve
eseguire prima un harness. Richiedono la build MPI seriale (un solo rank),
un thread OpenMP e l'integratore esplicito supportato.

Dalla radice del repository:

```sh
./ninja -C build SU2_CFD/src/SU2_CFD -j 4
cd TestCases/primitive_euler_examples/uniform
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD uniform.cfg
```

Dalla cartella `TestCases/primitive_euler_examples/vortex`:

```sh
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD vortex_regular.cfg
OMP_NUM_THREADS=1 ../../../build/SU2_CFD/src/SU2_CFD vortex_jitter.cfg
```

Il caso uniforme compie 8 passi da 0.01; i vortici N=16 compiono 512 passi
fino a T=0.5. I due vortici hanno prefissi di output distinti e possono
essere eseguiti nella stessa directory. Nessuna sorgente è attiva.

In ParaView aprire i file `uniform_*.vtm`, `vortex_regular_*.vtm` o
`vortex_jitter_*.vtm` come serie, premere Apply e selezionare il blocco
Internal. `Density`, `Velocity`, `Pressure`, `Physical_Time`, `Entropy`,
`Dual_Volume` e gli altri campi sono prodotti dal sistema nativo SU2.
Per confrontare pressioni mantenere la stessa scala cromatica.

Il writer XML nativo SU2 usa Float32; i file `*_state_*.dat` prodotti con
`OUTPUT_FILES=RESTART` conservano invece coordinate e W in Float64. Sono
checkpoint di output: il caricamento restart del primitivo non è supportato.
Il tempo e i residui di stadio sono nella history, con precisione 17.

Il bordo CUSTOM è uno **stato esterno** debole path-conservative, uniforme
o analitico. Non è una parete né una generica condizione caratteristica
di ingresso/uscita. Gli altri tipi di bordo non implementati sono rifiutati.

Le mesh sono riprese esattamente dagli archivi FASE 7; per generare altre
risoluzioni supportate (8,16,32,64) in una directory di lavoro:

```sh
python3 UnitTests/primitive_euler_driver_cases.py --output-dir /private/tmp/su2-phase8-examples --n 16
```

Configurazione, limiti, test e risultati sono nel
[rapporto FASE 8](../../Docs/primitive_euler_phase8_su2_integration.md).
