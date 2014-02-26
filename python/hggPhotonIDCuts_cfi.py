import FWCore.ParameterSet.Config as cms

hggPhotonIDCuts = cms.PSet(
#  ----------------------------sob value=0.0002     end of iteration 8          Loose  ---------------------------
      cutsubleadisosumoet0 = cms.vdouble(     8.2,       4.1,       5.4,       2.6),
   cutsubleadisosumoetbad0 = cms.vdouble(      67,        69,        85,       7.2),
    cutsubleadtrkisooetom0 = cms.vdouble(     7.5,       4.5,       5.2,       2.5),
          cutsubleadsieie0 = cms.vdouble(  0.0112,    0.0102,     0.029,     0.028),
         cutsubleadhovere0 = cms.vdouble(    0.09,     0.089,     0.101,     0.073),
             cutsubleadr90 = cms.vdouble(    0.94,      0.31,      0.92,      0.29),
  cutsublead_drtotk_25_990 = cms.vdouble(    0.26,     0.029,    0.0062,    0.0055),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.940332        fake=0.0848119        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0004     end of iteration 6          Medium  ---------------------------
      cutsubleadisosumoet1 = cms.vdouble(     6.4,       3.2,       3.4,       2.2),
   cutsubleadisosumoetbad1 = cms.vdouble(      64,      10.8,        13,       3.5),
    cutsubleadtrkisooetom1 = cms.vdouble(     6.4,       3.4,       3.8,       2.1),
          cutsubleadsieie1 = cms.vdouble(  0.0109,      0.01,     0.029,     0.028),
         cutsubleadhovere1 = cms.vdouble(   0.089,     0.079,      0.09,     0.061),
             cutsubleadr91 = cms.vdouble(    0.94,      0.32,      0.94,      0.29),
  cutsublead_drtotk_25_991 = cms.vdouble(    0.98,     0.029,    0.0109,    0.0111),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.918779        fake=0.0596949        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0008     end of iteration 6          Tight  ---------------------------
      cutsubleadisosumoet2 = cms.vdouble(     4.7,       2.8,       2.5,      1.46),
   cutsubleadisosumoetbad2 = cms.vdouble(      62,       5.2,       7.3,       2.5),
    cutsubleadtrkisooetom2 = cms.vdouble(     4.7,       2.9,       3.8,      1.63),
          cutsubleadsieie2 = cms.vdouble(  0.0107,    0.0099,     0.028,     0.027),
         cutsubleadhovere2 = cms.vdouble(   0.087,     0.065,     0.087,      0.05),
             cutsubleadr92 = cms.vdouble(    0.94,      0.34,      0.94,      0.29),
  cutsublead_drtotk_25_992 = cms.vdouble(       1,     0.029,     0.021,     0.028),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.887763        fake=0.042122        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0016     end of iteration 6          SuperTight  ---------------------------
      cutsubleadisosumoet3 = cms.vdouble(     3.8,       2.2,      1.77,      1.29),
   cutsubleadisosumoetbad3 = cms.vdouble(    11.7,       3.4,       3.9,      1.84),
    cutsubleadtrkisooetom3 = cms.vdouble(     3.5,       2.2,       2.3,      1.45),
          cutsubleadsieie3 = cms.vdouble(  0.0106,    0.0097,     0.028,     0.027),
         cutsubleadhovere3 = cms.vdouble(   0.082,     0.062,     0.065,     0.048),
             cutsubleadr93 = cms.vdouble(    0.94,      0.36,      0.94,      0.32),
  cutsublead_drtotk_25_993 = cms.vdouble(       1,     0.062,      0.97,      0.97),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.844291        fake=0.0290732        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0032     end of iteration 6          HyperTight1  ---------------------------
      cutsubleadisosumoet4 = cms.vdouble(     3.2,      1.76,      1.39,      1.18),
   cutsubleadisosumoetbad4 = cms.vdouble(     6.1,       2.7,       2.8,      0.66),
    cutsubleadtrkisooetom4 = cms.vdouble(     3.4,      1.86,      1.67,      1.44),
          cutsubleadsieie4 = cms.vdouble(  0.0104,    0.0094,     0.028,     0.025),
         cutsubleadhovere4 = cms.vdouble(   0.076,      0.03,     0.047,     0.046),
             cutsubleadr94 = cms.vdouble(    0.94,      0.41,      0.94,      0.34),
  cutsublead_drtotk_25_994 = cms.vdouble(       1,      0.97,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.774013        fake=0.0190779        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.00625     end of iteration 6          HyperTight2  ---------------------------
      cutsubleadisosumoet5 = cms.vdouble(     2.6,      1.31,      1.33,      0.82),
   cutsubleadisosumoetbad5 = cms.vdouble(     5.1,      1.62,      1.38,     -0.22),
    cutsubleadtrkisooetom5 = cms.vdouble(     2.9,       1.6,      1.55,      1.44),
          cutsubleadsieie5 = cms.vdouble(  0.0101,    0.0093,     0.027,     0.023),
         cutsubleadhovere5 = cms.vdouble(   0.048,    0.0189,     0.032,    0.0085),
             cutsubleadr95 = cms.vdouble(    0.94,      0.47,      0.94,      0.52),
  cutsublead_drtotk_25_995 = cms.vdouble(       1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.676391        fake=0.0121019        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0125     end of iteration 6          HyperTight3  ---------------------------
      cutsubleadisosumoet6 = cms.vdouble(    1.85,      0.96,      1.21,    -0.029),
   cutsubleadisosumoetbad6 = cms.vdouble(     3.7,      0.97,      1.38,     -0.88),
    cutsubleadtrkisooetom6 = cms.vdouble(    1.93,       1.4,      1.48,     0.056),
          cutsubleadsieie6 = cms.vdouble(  0.0099,    0.0092,     0.027,     0.023),
         cutsubleadhovere6 = cms.vdouble(   0.042,    0.0173,     0.023,    0.0085),
             cutsubleadr96 = cms.vdouble(    0.94,      0.69,      0.97,      0.52),
  cutsublead_drtotk_25_996 = cms.vdouble(       1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.524271        fake=0.00631764        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.025     end of iteration 6          HyperTight4  ---------------------------
      cutsubleadisosumoet7 = cms.vdouble(    1.31,       0.3,      1.15,    -0.029),
   cutsubleadisosumoetbad7 = cms.vdouble(    1.72,      0.69,      1.14,     -0.88),
    cutsubleadtrkisooetom7 = cms.vdouble(    1.42,      0.76,      1.48,     0.056),
          cutsubleadsieie7 = cms.vdouble(  0.0098,     0.009,     0.026,     0.023),
         cutsubleadhovere7 = cms.vdouble(   0.037,   0.00049,    0.0198,   0.00024),
             cutsubleadr97 = cms.vdouble(    0.94,      0.69,      0.97,      0.73),
  cutsublead_drtotk_25_997 = cms.vdouble(       1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.383546        fake=0.00362626        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.045     end of iteration 6          HyperTight5  ---------------------------
      cutsubleadisosumoet8 = cms.vdouble(    0.94,     0.112,     0.039,    -0.029),
   cutsubleadisosumoetbad8 = cms.vdouble(    0.86,      0.45,     -0.37,     -0.88),
    cutsubleadtrkisooetom8 = cms.vdouble(    1.21,      0.51,      0.27,     0.056),
          cutsubleadsieie8 = cms.vdouble(  0.0097,     0.009,     0.026,     0.023),
         cutsubleadhovere8 = cms.vdouble(   0.028,   1.4e-05,    0.0198,     7e-06),
             cutsubleadr98 = cms.vdouble(    0.94,      0.69,      0.97,      0.73),
  cutsublead_drtotk_25_998 = cms.vdouble(       1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.265878        fake=0.00218518        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.08     end of iteration 6          HyperTight6  ---------------------------
      cutsubleadisosumoet9 = cms.vdouble(     0.3,     0.072,     -0.04,     -0.97),
   cutsubleadisosumoetbad9 = cms.vdouble(    0.59,      0.42,     -0.84,     -1.77),
    cutsubleadtrkisooetom9 = cms.vdouble(    0.43,       0.4,    0.0075,     0.056),
          cutsubleadsieie9 = cms.vdouble(  0.0094,     0.009,     0.024,     0.023),
         cutsubleadhovere9 = cms.vdouble(  0.0071,         0,   0.00055,         0),
             cutsubleadr99 = cms.vdouble(    0.95,      0.69,      0.97,      0.84),
  cutsublead_drtotk_25_999 = cms.vdouble(       1,         1,         1,         1),

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.109842        fake=0.000923024        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.15     end of iteration 6          HyperTight7  ---------------------------
      cutsubleadisosumoet10 = cms.vdouble(   0.146,     0.058,     -0.04,     -1.16),
   cutsubleadisosumoetbad10 = cms.vdouble(    0.42,      0.41,     -0.84,     -1.77),
    cutsubleadtrkisooetom10 = cms.vdouble(    0.43,    0.0111,    0.0075,    0.0073),
          cutsubleadsieie10 = cms.vdouble(  0.0094,     0.009,     0.024,     0.023),
         cutsubleadhovere10 = cms.vdouble(  0.0002,         0,   1.6e-05,         0),
             cutsubleadr910 = cms.vdouble(    0.95,      0.69,      0.97,      0.85),
  cutsublead_drtotk_25_9910 = cms.vdouble(       1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.0898582        fake=0.00079285        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.15     end of iteration 6          Preselection  ---------------------------
      cutsubleadisosumoet11 = cms.vdouble(     7.0,       7.0,       4.0,       4.0),
   cutsubleadisosumoetbad11 = cms.vdouble(    100.,      100.,      100.,      100.),
    cutsubleadtrkisooetom11 = cms.vdouble(     7.0,       7.0,       4.0,       4.0),
          cutsubleadsieie11 = cms.vdouble(   0.012,     0.012,      0.03,      0.03),
         cutsubleadhovere11 = cms.vdouble(     99.,       99.,       99.,       99.),
             cutsubleadr911 = cms.vdouble(      0.,        0.,        0.,        0.),
  cutsublead_drtotk_25_9911 = cms.vdouble(       1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.0898582        fake=0.00079285        <<<<<<<<<<<<<<<<<<<<<<<<<<<<

###########################################################################
########################## 6 categories ###################################
###########################################################################
#   ----------------------------sob value=0.0002     end of iteration 8          Loose  ---------------------------
      cutsubleadisosumoet6c0 = cms.vdouble(    14.1,      11.7,       9.8,        11,       9.2,       8.9),
   cutsubleadisosumoetbad6c0 = cms.vdouble(      73,        59,        85,        94,        74,      14.4),
    cutsubleadtrkisooetom6c0 = cms.vdouble(     7.9,       5.6,       4.5,       5.9,       4.2,         3),
          cutsubleadsieie6c0 = cms.vdouble(  0.0114,     0.011,    0.0101,     0.029,     0.029,     0.028),
         cutsubleadhovere6c0 = cms.vdouble(   0.091,     0.079,     0.084,     0.103,      0.06,     0.099),
             cutsubleadr96c0 = cms.vdouble(    0.94,       0.9,      0.26,      0.94,       0.9,      0.28),
  cutsublead_drtotk_25_996c0 = cms.vdouble(    0.32,      0.99,    0.0097,    0.0064,    0.0048,    0.0061),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.947448        fake=0.0784026        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0004     end of iteration 6          Medium  ---------------------------
      cutsubleadisosumoet6c1 = cms.vdouble(    12.4,      10.2,       9.2,        10,       8.3,       8.7),
   cutsubleadisosumoetbad6c1 = cms.vdouble(      71,        50,        22,        24,      18.8,      12.7),
    cutsubleadtrkisooetom6c1 = cms.vdouble(     6.6,       4.7,       4.4,       5.7,       2.4,       2.2),
          cutsubleadsieie6c1 = cms.vdouble(   0.011,    0.0107,    0.0098,     0.029,     0.029,     0.027),
         cutsubleadhovere6c1 = cms.vdouble(   0.089,     0.076,      0.08,     0.091,     0.051,     0.063),
             cutsubleadr96c1 = cms.vdouble(    0.94,       0.9,      0.27,      0.94,       0.9,      0.28),
  cutsublead_drtotk_25_996c1 = cms.vdouble(    0.98,         1,     0.012,    0.0111,    0.0086,    0.0122),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.928068        fake=0.0572766        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0008     end of iteration 6          Tight  ---------------------------
      cutsubleadisosumoet6c2 = cms.vdouble(    11.2,       9.3,         9,       9.2,       7.9,       8.2),
   cutsubleadisosumoetbad6c2 = cms.vdouble(      71,      16.8,      13.7,      15.6,      13.2,      10.4),
    cutsubleadtrkisooetom6c2 = cms.vdouble(     5.8,         4,       2.9,       4.6,       2.3,      1.93),
          cutsubleadsieie6c2 = cms.vdouble(  0.0108,    0.0105,    0.0097,     0.028,     0.028,     0.027),
         cutsubleadhovere6c2 = cms.vdouble(   0.083,      0.07,     0.061,     0.088,      0.05,     0.061),
             cutsubleadr96c2 = cms.vdouble(    0.94,       0.9,      0.32,      0.94,       0.9,      0.28),
  cutsublead_drtotk_25_996c2 = cms.vdouble(       1,         1,    0.0146,     0.021,      0.36,     0.026),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.895288        fake=0.0392946        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0016     end of iteration 6          SuperTight  ---------------------------
      cutsubleadisosumoet6c3 = cms.vdouble(      10,       8.8,       8.7,       8.5,       7.9,       7.5),
   cutsubleadisosumoetbad6c3 = cms.vdouble(      21,      14.4,      11.5,      12.9,      10.2,       9.7),
    cutsubleadtrkisooetom6c3 = cms.vdouble(       4,       3.4,       2.5,      1.89,       1.7,      1.67),
          cutsubleadsieie6c3 = cms.vdouble(  0.0106,    0.0102,    0.0095,     0.028,     0.028,     0.026),
         cutsubleadhovere6c3 = cms.vdouble(   0.082,      0.06,     0.049,     0.077,      0.05,     0.061),
             cutsubleadr96c3 = cms.vdouble(    0.94,       0.9,      0.32,      0.94,       0.9,      0.28),
  cutsublead_drtotk_25_996c3 = cms.vdouble(       1,         1,     0.022,      0.97,      0.98,      0.97),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.849144        fake=0.0256315        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0032     end of iteration 6          HyperTight1  ---------------------------
      cutsubleadisosumoet6c4 = cms.vdouble(     9.3,       8.1,       7.5,       7.9,       6.3,       6.6),
   cutsubleadisosumoetbad6c4 = cms.vdouble(    15.4,      10.8,      10.2,      11.5,       8.5,       8.8),
    cutsubleadtrkisooetom6c4 = cms.vdouble(     3.6,       2.9,      1.96,      1.63,      1.43,      1.36),
          cutsubleadsieie6c4 = cms.vdouble(  0.0105,      0.01,    0.0094,     0.028,     0.027,     0.026),
         cutsubleadhovere6c4 = cms.vdouble(   0.068,      0.06,     0.027,     0.055,     0.047,     0.047),
             cutsubleadr96c4 = cms.vdouble(    0.94,       0.9,      0.35,      0.94,       0.9,       0.3),
  cutsublead_drtotk_25_996c4 = cms.vdouble(       1,         1,     0.091,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.784975        fake=0.0163432        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.00625     end of iteration 6          HyperTight2  ---------------------------
      cutsubleadisosumoet6c5 = cms.vdouble(     8.6,       6.6,       6.9,       7.1,         6,       6.3),
   cutsubleadisosumoetbad6c5 = cms.vdouble(    13.4,      10.4,       9.2,       9.3,         7,       7.3),
    cutsubleadtrkisooetom6c5 = cms.vdouble(     2.7,       2.9,      1.56,      1.56,      1.06,      1.36),
          cutsubleadsieie6c5 = cms.vdouble(  0.0104,    0.0098,    0.0094,     0.028,     0.027,     0.024),
         cutsubleadhovere6c5 = cms.vdouble(   0.058,     0.023,    0.0174,     0.055,     0.047,     0.047),
             cutsubleadr96c5= cms.vdouble(    0.94,       0.9,      0.37,      0.94,       0.9,      0.58),
  cutsublead_drtotk_25_996c5 = cms.vdouble(       1,         1,      0.97,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.69497        fake=0.00960486        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0125     end of iteration 6          HyperTight3  ---------------------------
      cutsubleadisosumoet6c6 = cms.vdouble(       8,       6.6,       6.5,       5.3,       4.9,       5.7),
   cutsubleadisosumoetbad6c6 = cms.vdouble(    11.4,       8.9,       8.3,       8.1,         7,       6.7),
    cutsubleadtrkisooetom6c6 = cms.vdouble(     2.3,       2.6,      1.54,      1.56,      1.06,      1.07),
          cutsubleadsieie6c6 = cms.vdouble(  0.0101,    0.0097,    0.0093,     0.028,     0.025,     0.023),
         cutsubleadhovere6c6 = cms.vdouble(   0.038,     0.022,     0.014,     0.043,    0.0116,     0.047),
             cutsubleadr96c6 = cms.vdouble(    0.94,       0.9,      0.37,      0.95,       0.9,      0.58),
  cutsublead_drtotk_25_996c6 = cms.vdouble(       1,         1,         1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.570812        fake=0.005191        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.025     end of iteration 6          HyperTight4  ---------------------------
      cutsubleadisosumoet6c7 = cms.vdouble(     6.8,       6.1,         5,       4.8,       4.6,       4.9),
   cutsubleadisosumoetbad6c7 = cms.vdouble(     9.2,       8.1,       7.7,       7.4,         7,       6.7),
    cutsubleadtrkisooetom6c7 = cms.vdouble(     2.1,      1.79,      1.07,      1.56,       0.8,      0.26),
          cutsubleadsieie6c7 = cms.vdouble(  0.0099,    0.0092,    0.0092,     0.028,     0.025,     0.022),
         cutsubleadhovere6c7 = cms.vdouble(   0.032,    0.0154,    0.0111,     0.023,    0.0107,    0.0138),
             cutsubleadr96c7 = cms.vdouble(    0.94,       0.9,      0.39,      0.95,       0.9,      0.69),
  cutsublead_drtotk_25_996c7 = cms.vdouble(       1,         1,         1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.415682        fake=0.0022302        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.045     end of iteration 6          HyperTight5  ---------------------------
      cutsubleadisosumoet6c8 = cms.vdouble(     6.5,       5.7,       4.8,       4.4,       4.6,       4.5),
   cutsubleadisosumoetbad6c8 = cms.vdouble(     8.3,       7.7,       6.4,       6.8,       5.3,       6.5),
    cutsubleadtrkisooetom6c8 = cms.vdouble(     2.1,      1.55,      0.25,       0.8,      0.48,      0.26),
          cutsubleadsieie6c8 = cms.vdouble(  0.0098,    0.0092,    0.0088,     0.025,     0.025,     0.022),
         cutsubleadhovere6c8 = cms.vdouble(   0.026,    0.0137,    0.0111,    0.0192,    0.0073,    0.0083),
             cutsubleadr96c8 = cms.vdouble(    0.94,       0.9,      0.52,      0.95,      0.92,      0.76),
  cutsublead_drtotk_25_996c8 = cms.vdouble(       1,         1,         1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.299123        fake=0.00101069        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.08     end of iteration 6          HyperTight6  ---------------------------
      cutsubleadisosumoet6c9 = cms.vdouble(     6.1,       5.6,       4.8,       4.4,       4.6,       4.5),
   cutsubleadisosumoetbad6c9 = cms.vdouble(     7.8,       7.4,       6.2,       5.4,       5.3,       6.5),
    cutsubleadtrkisooetom6c9 = cms.vdouble(    1.83,       1.2,     0.029,      0.68,      0.48,      0.26),
          cutsubleadsieie6c9 = cms.vdouble(  0.0098,    0.0091,    0.0086,     0.024,     0.023,     0.022),
         cutsubleadhovere6c9 = cms.vdouble(   0.026,    0.0137,   0.00114,   0.00068,    0.0073,    0.0083),
             cutsubleadr96c9 = cms.vdouble(    0.94,       0.9,      0.52,      0.95,      0.92,      0.76),
  cutsublead_drtotk_25_996c9 = cms.vdouble(       1,         1,         1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.227233        fake=0.000514037        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.15     end of iteration 6          HyperTight7  ---------------------------
      cutsubleadisosumoet6c10= cms.vdouble(     5.2,       5.6,       4.8,       4.4,       4.6,       4.5),
   cutsubleadisosumoetbad6c10= cms.vdouble(     7.7,       7.1,       6.1,       5.3,       5.2,       6.5),
    cutsubleadtrkisooetom6c10= cms.vdouble(     1.6,     0.033,   0.00082,      0.68,      0.48,      0.26),
          cutsubleadsieie6c10= cms.vdouble(  0.0098,    0.0091,    0.0086,     0.024,     0.023,     0.022),
         cutsubleadhovere6c10= cms.vdouble(   0.026,    0.0082,   3.2e-05,   3.2e-05,    0.0073,    0.0083),
             cutsubleadr96c10= cms.vdouble(    0.94,       0.9,      0.52,      0.95,      0.92,      0.76),
  cutsublead_drtotk_25_996c10= cms.vdouble(       1,         1,         1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.167242        fake=0.000199112        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.15     end of iteration 6          Preselection ---------------------------
      cutsubleadisosumoet6c11= cms.vdouble(     7.0,       7.0,       7.0,       4.0,       4.0,       4.0),
   cutsubleadisosumoetbad6c11= cms.vdouble(    100.,      100.,      100.,      100.,      100.,      100.),
    cutsubleadtrkisooetom6c11= cms.vdouble(     7.0,       7.0,       7.0,       4.0,       4.0,       4.0),
          cutsubleadsieie6c11= cms.vdouble(   0.012,     0.012,     0.012,      0.03,      0.03,      0.03),
         cutsubleadhovere6c11= cms.vdouble(     99.,       99.,       99.,       99.,       99.,       99.),
             cutsubleadr96c11= cms.vdouble(      0.,        0.,        0.,        0.,        0.,        0.),
  cutsublead_drtotk_25_996c11= cms.vdouble(       1,         1,         1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.167242        fake=0.000199112        <<<<<<<<<<<<<<<<<<<<<<<<<<<<

###########################################################################
########################## 6 categories PF ################################
###########################################################################
#  ----------------------------sob value=0.000222222     end of iteration 9          Loose  ---------------------------
    cutsubleadpfisosumoet0 = cms.vdouble(     9.2,       7.4,         6,       9.3,        10,         7),
 cutsubleadpfisosumoetbad0 = cms.vdouble(      47,        55,      11.3,      19.3,        58,       7.7),
       cutsubleadpfchisooet0 = cms.vdouble(     6.3,       4.8,       4.4,       4.4,       7.1,       4.2),
#   cutsubleadpfhcalisooetom0 = cms.vdouble(    11.6,        22,       8.7,       8.3,       8.9,       7.5),
          cutsubleadpfsieie0 = cms.vdouble(  0.0116,    0.0114,    0.0103,     0.032,     0.031,      0.03),
         cutsubleadpfhovere0 = cms.vdouble(   0.129,     0.143,     0.112,     0.146,     0.136,     0.147),
             cutsubleadpfr90 = cms.vdouble(    0.94,       0.9,      0.23,      0.94,       0.9,      0.23),
  cutsubleadpf_drtotk_25_990 = cms.vdouble(    0.99,      0.96,    0.0064,    0.0066,    0.0041,    0.0032),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.969349        fake=0.702721        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.000388889     end of iteration 8          Medium  ---------------------------
    cutsubleadpfisosumoet1 = cms.vdouble(       8,       6.5,       5.1,       6.2,       6.7,       4.8),
 cutsubleadpfisosumoetbad1 = cms.vdouble(      44,      11.4,       8.2,      10.1,       9.5,       6.5),
       cutsubleadpfchisooet1 = cms.vdouble(     5.6,       3.7,       3.5,       3.7,       4.5,       3.4),
#   cutsubleadpfhcalisooetom1 = cms.vdouble(    11.6,        22,       7.8,       8.2,       6.5,       7.4),
          cutsubleadpfsieie1 = cms.vdouble(  0.0116,    0.0114,    0.0101,     0.029,      0.03,     0.029),
         cutsubleadpfhovere1 = cms.vdouble(   0.129,     0.142,      0.11,     0.146,     0.136,     0.075),
             cutsubleadpfr91 = cms.vdouble(    0.94,       0.9,      0.23,      0.94,       0.9,      0.23),
  cutsubleadpf_drtotk_25_991 = cms.vdouble(       1,         1,    0.0114,    0.0081,    0.0068,    0.0109),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.949368        fake=0.617595        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.000555556     end of iteration 8          Tight  ---------------------------
    cutsubleadpfisosumoet2 = cms.vdouble(     7.2,       6.2,       4.7,       5.7,       6.7,       4.3),
 cutsubleadpfisosumoetbad2 = cms.vdouble(    11.8,       8.2,       6.8,       8.4,       6.5,       5.6),
       cutsubleadpfchisooet2 = cms.vdouble(     4.7,       3.7,       2.5,       3.4,       4.5,       2.3),
#   cutsubleadpfhcalisooetom2 = cms.vdouble(    11.6,        22,       7.8,       8.2,       6.4,       7.4),
          cutsubleadpfsieie2 = cms.vdouble(  0.0113,    0.0111,    0.0101,     0.028,      0.03,     0.028),
         cutsubleadpfhovere2 = cms.vdouble(   0.123,     0.142,      0.11,     0.146,     0.134,     0.061),
             cutsubleadpfr92 = cms.vdouble(    0.94,       0.9,      0.24,      0.94,       0.9,      0.23),
  cutsubleadpf_drtotk_25_992 = cms.vdouble(       1,         1,    0.0125,    0.0094,      0.02,    0.0145),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.928023        fake=0.559129        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.000722222     end of iteration 8          SuperTight  ---------------------------
    cutsubleadpfisosumoet3 = cms.vdouble(     6.1,       5.5,       4.7,       5.6,       6.3,       3.8),
 cutsubleadpfisosumoetbad3 = cms.vdouble(    10.9,       7.5,       6.2,       6.9,       5.6,       5.1),
       cutsubleadpfchisooet3 = cms.vdouble(       4,       3.3,       2.4,       3.2,       4.3,         2),
#   cutsubleadpfhcalisooetom3 = cms.vdouble(    11.5,        22,       7.8,         8,       6.2,       7.4),
          cutsubleadpfsieie3 = cms.vdouble(  0.0112,    0.0108,      0.01,     0.028,     0.029,     0.027),
         cutsubleadpfhovere3 = cms.vdouble(   0.123,     0.142,      0.11,     0.146,     0.108,     0.059),
             cutsubleadpfr93 = cms.vdouble(    0.94,       0.9,      0.24,      0.94,       0.9,      0.23),
  cutsubleadpf_drtotk_25_993 = cms.vdouble(       1,         1,      0.06,    0.0163,     0.057,    0.0168),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.903277        fake=0.511431        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.000944444     end of iteration 8          HyperTight1  ---------------------------
    cutsubleadpfisosumoet4 = cms.vdouble(     5.6,       4.8,       4.3,       5.6,       3.9,         3),
 cutsubleadpfisosumoetbad4 = cms.vdouble(     9.8,       6.7,       5.6,       5.3,       4.6,         4),
       cutsubleadpfchisooet4 = cms.vdouble(     3.4,       2.3,       2.4,       2.4,       1.9,      1.35),
#   cutsubleadpfhcalisooetom4 = cms.vdouble(    11.5,        22,       7.8,       4.3,       5.6,       3.3),
          cutsubleadpfsieie4 = cms.vdouble(  0.0108,    0.0107,    0.0099,     0.028,     0.028,     0.027),
         cutsubleadpfhovere4 = cms.vdouble(   0.123,     0.142,      0.11,     0.146,      0.09,     0.059),
             cutsubleadpfr94 = cms.vdouble(    0.94,       0.9,      0.38,      0.94,       0.9,      0.23),
  cutsubleadpf_drtotk_25_994 = cms.vdouble(       1,         1,      0.44,      0.36,      0.39,    0.0182),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.870999        fake=0.459021        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.00122222     end of iteration 8          HyperTight2  ---------------------------
    cutsubleadpfisosumoet5 = cms.vdouble(     5.2,       4.6,       3.6,         4,       3.2,       2.5),
 cutsubleadpfisosumoetbad5 = cms.vdouble(       9,       6.4,       4.6,       4.9,       4.2,       3.8),
       cutsubleadpfchisooet5 = cms.vdouble(     3.2,       2.1,       2.1,      1.97,      1.33,     0.028),
#   cutsubleadpfhcalisooetom5 = cms.vdouble(    11.5,        21,       7.8,       4.2,      1.93,       2.9),
          cutsubleadpfsieie5 = cms.vdouble(  0.0107,    0.0107,    0.0098,     0.028,     0.027,     0.027),
         cutsubleadpfhovere5 = cms.vdouble(   0.123,     0.142,      0.11,     0.108,     0.068,     0.059),
             cutsubleadpfr95 = cms.vdouble(    0.94,       0.9,      0.38,      0.94,       0.9,      0.23),
  cutsubleadpf_drtotk_25_995 = cms.vdouble(       1,         1,      0.99,      0.99,      0.96,     0.019),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.830285        fake=0.404951        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.00166667     end of iteration 8          HyperTight3  ---------------------------
    cutsubleadpfisosumoet6 = cms.vdouble(     4.8,       4.3,       2.7,       2.7,       2.5,       2.4),
 cutsubleadpfisosumoetbad6 = cms.vdouble(     6.2,         6,       4.2,         4,       4.2,       3.4),
       cutsubleadpfchisooet6 = cms.vdouble(     2.7,       2.1,      1.68,      1.25,       1.1,     0.028),
#   cutsubleadpfhcalisooetom6 = cms.vdouble(    11.5,        21,       7.6,       4.2,      1.93,      1.83),
          cutsubleadpfsieie6 = cms.vdouble(  0.0107,    0.0105,    0.0098,     0.027,     0.027,     0.027),
         cutsubleadpfhovere6 = cms.vdouble(   0.123,     0.122,      0.11,     0.108,     0.026,     0.059),
             cutsubleadpfr96 = cms.vdouble(    0.94,       0.9,      0.38,      0.94,       0.9,      0.23),
  cutsubleadpf_drtotk_25_996 = cms.vdouble(       1,         1,         1,         1,         1,      0.53),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.727129        fake=0.305297        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.00222222     end of iteration 8          HyperTight4  ---------------------------
    cutsubleadpfisosumoet7 = cms.vdouble(     2.8,       2.8,       2.4,       2.5,       2.4,       2.3),
 cutsubleadpfisosumoetbad7 = cms.vdouble(     6.2,       5.4,       3.5,         4,       4.2,       2.6),
       cutsubleadpfchisooet7 = cms.vdouble(    1.32,       1.3,      1.45,      0.76,      0.97,   0.00047),
#   cutsubleadpfhcalisooetom7 = cms.vdouble(     9.5,       4.3,       2.4,       4.2,      1.93,      1.34),
          cutsubleadpfsieie7 = cms.vdouble(  0.0107,    0.0105,    0.0098,     0.027,     0.027,     0.027),
         cutsubleadpfhovere7 = cms.vdouble(    0.09,     0.091,      0.11,     0.067,    0.0198,     0.053),
             cutsubleadpfr97 = cms.vdouble(    0.94,       0.9,      0.38,      0.94,       0.9,      0.32),
  cutsubleadpf_drtotk_25_997 = cms.vdouble(       1,         1,         1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.534842        fake=0.174365        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.00305556     end of iteration 8          HyperTight5  ---------------------------
    cutsubleadpfisosumoet8 = cms.vdouble(     2.5,       2.4,       2.2,       2.3,       2.2,      1.98),
 cutsubleadpfisosumoetbad8 = cms.vdouble(     5.7,       4.7,       2.7,         4,       4.2,       2.1),
       cutsubleadpfchisooet8 = cms.vdouble(    1.07,         1,      1.34,      0.69,      0.69,     5e-06),
#   cutsubleadpfhcalisooetom8 = cms.vdouble(     9.5,       4.3,       2.4,       2.4,      1.93,      1.17),
          cutsubleadpfsieie8 = cms.vdouble(  0.0107,    0.0104,    0.0097,     0.027,     0.024,     0.025),
         cutsubleadpfhovere8 = cms.vdouble(   0.086,     0.091,     0.059,     0.057,    0.0198,     0.036),
             cutsubleadpfr98 = cms.vdouble(    0.94,       0.9,      0.42,      0.95,       0.9,      0.53),
  cutsubleadpf_drtotk_25_998 = cms.vdouble(       1,         1,         1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.405524        fake=0.110157        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.00416667     end of iteration 8          HyperTight6  ---------------------------
    cutsubleadpfisosumoet9 = cms.vdouble(     2.4,       2.2,      1.88,       2.2,       2.2,       1.7),
 cutsubleadpfisosumoetbad9 = cms.vdouble(     5.6,       3.6,       2.1,       2.4,       4.2,      1.55),
       cutsubleadpfchisooet9 = cms.vdouble(    0.98,      0.68,      1.34,     0.007,      0.52,         0),
#   cutsubleadpfhcalisooetom9 = cms.vdouble(     5.6,       2.5,      1.27,         2,      1.02,      1.17),
          cutsubleadpfsieie9 = cms.vdouble(  0.0107,    0.0103,    0.0097,     0.025,     0.023,     0.024),
         cutsubleadpfhovere9 = cms.vdouble(   0.039,      0.04,     0.059,     0.032,    0.0198,    0.0178),
             cutsubleadpfr99 = cms.vdouble(    0.94,       0.9,      0.42,      0.95,       0.9,      0.67),
  cutsubleadpf_drtotk_25_999 = cms.vdouble(       1,         1,         1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.288559        fake=0.0665108        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.00555556     end of iteration 8          HyperTight7  ---------------------------
   cutsubleadpfisosumoet10 = cms.vdouble(     2.2,       2.1,      1.37,       2.1,       2.2,       1.7),
cutsubleadpfisosumoetbad10 = cms.vdouble(     4.9,       3.6,      1.17,       2.1,       3.1,      1.55),
      cutsubleadpfchisooet10 = cms.vdouble(    0.96,     0.149,    0.0135,     7e-05,     0.159,         0),
#  cutsubleadpfhcalisooetom10 = cms.vdouble(     2.1,      1.12,      0.69,      1.34,      0.78,       0.7),
         cutsubleadpfsieie10 = cms.vdouble(  0.0107,    0.0101,     0.009,     0.024,     0.023,     0.023),
        cutsubleadpfhovere10 = cms.vdouble(    0.02,    0.0164,     0.021,     0.029,     0.018,    0.0178),
            cutsubleadpfr910 = cms.vdouble(    0.94,       0.9,      0.46,      0.95,       0.9,      0.76),
 cutsubleadpf_drtotk_25_9910 = cms.vdouble(       1,         1,         1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.193294        fake=0.0391879        <<<<<<<<<<<<<<<<<<<<<<<<<<<<
#  ----------------------------sob value=0.0075     end of iteration 8          Preselection ---------------------------
      cutsubleadpfisosumoet11 = cms.vdouble(      8.,        8.,        8.,        7.,        7.,        7.),
   cutsubleadpfisosumoetbad11 = cms.vdouble(    100.,      100.,      100.,      100.,      100.,      100.),
       cutsubleadpfchisooet11 = cms.vdouble(      6.,        6.,        6.,        5.,        5.,        5.),
          cutsubleadpfsieie11 = cms.vdouble(   0.012,     0.012,     0.012,      0.03,      0.03,      0.03),
         cutsubleadpfhovere11 = cms.vdouble(     99.,       99.,       99.,       99.,       99.,       99.),
             cutsubleadpfr911 = cms.vdouble(      0.,        0.,        0.,        0.,        0.,        0.),
  cutsubleadpf_drtotk_25_9911 = cms.vdouble(       1,         1,         1,         1,         1,         1),
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  eff=0.0573526        fake=0.0108523        <<<<<<<<<<<<<<<<<<<<<<<<<<<<

)
