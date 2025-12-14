


              SIMIND Monte Carlo Simulation Program    V8.0  
------------------------------------------------------------------------------
 Phantom S : h2o       Crystal...: nai       InputFile.: spect_simulation  
 Phantom B : h2o       BackScatt.: pmt       OutputFile: spect_simulation_k
 Collimator: pb_sb2    SourceRout: smap      SourceImg.: pbpk_kidney_act_av
 Cover.....: al        ScoreRout.: scattwin  DensityImg: spect_preprocessin
------------------------------------------------------------------------------
 PhotonEnergy.......: 200          lu177     PhotonsPerProj....: 2399           
 EnergyResolution...: 8.8          Spectra   Activity..........: 1              
 MaxScatterOrder....: 3            si-me     DetectorLenght....: 21.82          
 DetectorWidth......: 53.3         SPECT     DetectorHeight....: 0.953          
 UpperEneWindowTresh: 220          x-rays    Distance to det...: 35.203         
 LowerEneWindowTresh: 180          BScatt    ShiftSource X.....: 0              
 PixelSize  I.......: 0.5          Random    ShiftSource Y.....: 0              
 PixelSize  J.......: 0.5          Cover     ShiftSource Z.....: 0              
 HalfLength S.......: 10.90908     Phantom   HalfLength P......: 10.90908       
 HalfWidth  S.......: 16           Resolut   HalfWidth  P......: 16             
 HalfHeight S.......: 16           Header    HalfHeight P......: 16             
 SourceType.........: F4 XcatMap   SaveMap   PhantomType.......: F4 XcatMap   
------------------------------------------------------------------------------
 GENERAL DATA
 keV/channel........: 0.5                    CutoffEnergy......: 0              
 Photons/Bq.........: 0.2264                 StartingAngle.....: 0              
 CameraOffset X.....: 0                      CoverThickness....: 0.1            
 CameraOffset Y.....: 0                      BackscatterThickn.: 10             
 MatrixSize I.......: 128                    IntrinsicResolut..: 0.34           
 MatrixSize J.......: 43                     AcceptanceAngle...: 4.13771        
 Emission type......: 2                      Initial Weight....: 94.37182       
 NN ScalingFactor...: 1.02768                Energy Channels...: 512            
                                                                              
 SPECT DATA
 RotationMode.......: -360                   Nr of Projections.: 64             
 RotationAngle......: 5.625                  Projection.[start]: 1              
 Orbital fraction...: 0                      Projection...[end]: 64             
 Center of Rotation File: spect_simulation_kidney_2.cor
                                                                              
 COLLIMATOR DATA FOR ROUTINE: Analytical          
 CollimatorCode.....: si-me                  CollimatorType....: Parallel 
 HoleSize X.........: 0.294                  Distance X........: 0.114          
 HoleSize Y.........: 0.33948                Distance Y........: 0.26847        
 CenterShift X......: 0.204                  X-Ray flag........: T              
 CenterShift Y......: 0.35334                CollimThickness...: 4.064          
 HoleShape..........: Hexagonal              Space Coll2Det....: 0              
 CollDepValue [57]..: 0                      CollDepValue [58].: 0              
 CollDepValue [59]..: 0                      CollDepValue [60].: 0              
                                                                              
 IMAGE-BASED PHANTOM DATA
 RotationCentre.....: 129,129                Bone definition...: 1170           
 CT-Pixel size......: 0.23976                Slice thickness...: 0.23976        
 StartImage.........: 1                      No of CT-Images...: 91             
 MatrixSize I.......: 256                    CTmapOrientation..: 0              
 MatrixSize J.......: 256                    StepSize..........: 0.1            
 CenterPoint I......: 129                    ShiftPhantom X....: 0              
 CenterPoint J......: 129                    ShiftPhantom Y....: 0              
 CenterPoint K......: 46.5                   ShiftPhantom Z....: 0              
                                                                              
 INFO FOR TCT file
 MatrixSize I.......: 128                    MatrixSize J......: 128            
 MatrixSize K.......: 43                     Units.............: mu                  
 Scout File.........: F
                                                                              
------------------------------------------------------------------------------
  Scattwin results: Window file: spect_simulation.win
  
  Win  WinAdded  Range(keV)   ScaleFactor
   1       0    169.4 - 187.2   1.000
   2       0    187.2 - 228.8   1.000
   3       0    228.8 - 252.9   1.000
  
  Win    Total    Scatter   Primary  S/P-Ratio S/T Ratio  Cps/MBq
   1   0.468E+02 0.443E+02 0.247E+01 0.180E+02 0.947E+00 0.731E+00
   2   0.144E+03 0.515E+02 0.925E+02 0.556E+00 0.357E+00 0.225E+01
   3   0.187E+01 0.790E+00 0.108E+01 0.733E+00 0.423E+00 0.292E-01
  
  Win  Geo(Air)  Pen(Air)  Sca(Air)  Geo(Tot)  Pen(Tot)  Sca(Tot)
   1   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
   2   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
   3   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
  
  Win   SC 1  SC 2  SC 3
   1   62.1% 31.2%  6.7%
   2   84.9% 13.3%  1.8%
   3   75.5% 14.1% 10.4%
                                                                              
 INTERACTIONS IN THE CRYSTAL
 MaxValue spectrum..: 6.450          
 MaxValue projection: 0.3580         
 CountRate spectrum.: 14.86          
 CountRate E-Window.: 2.483          
                                                                              
 SCATTER IN ENERGY WINDOW
 Scatter/Primary....: 0.76142        
 Scatter/Total......: 0.43228        
 Scatter order 1....: 80.13 %        
 Scatter order 2....: 16.98 %        
 Scatter order 3....: 2.89 %         
                                                                              
 CALCULATED DETECTOR PARAMETERS
 Efficiency E-window: 0.1513         
 Efficiency spectrum: 0.9055         
 Sensitivity Cps/MBq: 2.4831         
 Sensitivity Cpm/uCi: 5.5125         
                                                                              
 Simulation started.: 2025:12:13 00:41:09
 Simulation stopped.: 2025:12:13 00:41:17
 Elapsed time.......: 0 h, 0 m and 8 s
 DetectorHits.......: 6936           
 DetectorHits/CPUsec: 1097           
                                                                              
 OTHER INFORMATION
 XCAT Simulation
 Compiled 2025:10:06 with intel Linux 
 Current random number generator: ranmar
 Energy resolution as function of 1/sqrt(E)
 Header file: spect_simulation_kidney_2.h00
 Linear angle sampling within acceptance angle
 Inifile: simind.ini
 Command: spect_simulation spect_simulation_kidney_2 /fd:spect_preprocessing_atn_av.bin/fs:pbpk_kidney_act_av.bin/in:x22,3x/nn:1.027678370475769/cc:si-me/fi:lu177/02:10.909079706668853/05:10.909079706668853/08:21.82/10:53.30/14:-7/15:-7/28:0.5/29:64/31:0.23976001739501954/34:91/42:-10/76:128/77:43.636318826675414/78:256/79:256/rr:2
