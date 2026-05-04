# print.quaqc output matches snapshot

    Code
      print(FIXTURE)
    Output
      quaqc v1.1
      Run title: ---
      Run date: 2024-06-13 12:13:57 CEST
      Number of samples: 1
      Number of fails:   0
      Reports:
          [SUCCESS] SRR26098097.Chr4MtPt.bam
      To examine run data: $metadata
      To examine reports:  $reports
      See ?melt_reports or ?parse_quaqc for ways to access data.

# print.quaqc_report output matches snapshot

    Code
      print(FIXTURE$reports[[1]])
    Output
      Sample: SRR26098097.Chr4MtPt.bam
      Status: success
      Processed 1 sequences (18293156 bp).
      Reads: 9846473 total; 3867479 passing.
      GC histo: yes. Depths histo: yes. 
      Peak stats: yes. TSS pileup: yes.
      Available slots:
      $sample      -- Sample filename
      $success     -- Run status
      $params      -- Run parameters
      $genome      -- Genome stats
      $unfiltered  -- Pre-filter read stats
      $filtered    -- Post-filter read stats

