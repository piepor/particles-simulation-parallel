**Algoritmi paralleli provati**

*  ForceCompt_par (parallelo) - ComptPopulation_par (parallelo) - DumpPoplation (seriale, teoricamente eseguita in parallelo dalla CPU) - PopulationScreen (seriale, idem)
* ForceCompt_par (parallelo) - ReduceForce (parallelo) - ForceCompt_par (parallelo) - ReduceForce (parallelo) - Compt Population_par (parallelo) - DumpPopulation (seriale, teroicamente in parallelo dalla CPU) - ParticleScreen (seriale, idem)
* ForceCompt_par (parallelo) - ReduceForce (parallelo) - ForceCompt_par (parallelo) - ReduceForce (parallelo) - ForceCompt_par (parallelo) - ReduceForce (parallelo) - ForceCompt_par (parallelo) - ReduceForce (parallelo) - Compt Population_par (parallelo) - DumpPopulation (seriale, teroicamente in parallelo dalla CPU) - ParticleScreen (seriale, idem)
* ForceCompt_par (parallelo) - ComptPopulation_par (parallelo) - ParticleScreen_par (parallelo) - DumpPopulation (seriale, teroicamente in parallelo dalla CPU) - DumpGridValues (seriale, teroicamente in parallelo dalla CPU) - WriteImage (seriale, dopo tutto)



**Cose strane**

* latency report1.qdrep 8.574 micro secondi con occupancy teorica 50% mentre latency report5.qdrep 41.43 microsecondi con occupancy teorica 100%
* ReduceForce report5.qdrep latency 360 milli secondi con occpuancy teorica 50%
* memcpy si sovrappone? sembra di si sulla riga CUDA API invece sembra di no sotto CUDA HW

