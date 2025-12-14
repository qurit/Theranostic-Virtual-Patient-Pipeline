


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
 Center of Rotation File: spect_simulation_kidney_0.cor
                                                                              
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
   1   0.463E+02 0.439E+02 0.239E+01 0.184E+02 0.948E+00 0.723E+00
   2   0.143E+03 0.518E+02 0.915E+02 0.566E+00 0.361E+00 0.224E+01
   3   0.219E+01 0.110E+01 0.110E+01 0.100E+01 0.500E+00 0.343E-01
  
  Win  Geo(Air)  Pen(Air)  Sca(Air)  Geo(Tot)  Pen(Tot)  Sca(Tot)
   1   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
   2   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
   3   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
  
  Win   SC 1  SC 2  SC 3
   1   63.7% 30.5%  5.8%
   2   84.5% 13.8%  1.7%
   3   71.3% 23.0%  5.8%
                                                                              
 INTERACTIONS IN THE CRYSTAL
 MaxValue spectrum..: 6.680          
 MaxValue projection: 0.4826E-01     
 CountRate spectrum.: 14.87          
 CountRate E-Window.: 2.462          
                                                                              
 SCATTER IN ENERGY WINDOW
 Scatter/Primary....: 0.768          
 Scatter/Total......: 0.43439        
 Scatter order 1....: 80.82 %        
 Scatter order 2....: 16.93 %        
 Scatter order 3....: 2.25 %         
                                                                              
 CALCULATED DETECTOR PARAMETERS
 Efficiency E-window: 0.1496         
 Efficiency spectrum: 0.9034         
 Sensitivity Cps/MBq: 2.4618         
 Sensitivity Cpm/uCi: 5.4651         
                                                                              
 Simulation started.: 2025:12:13 00:41:09
 Simulation stopped.: 2025:12:13 00:41:18
 Elapsed time.......: 0 h, 0 m and 9 s
 DetectorHits.......: 7004           
 DetectorHits/CPUsec: 1084           
                                                                              
 OTHER INFORMATION
 XCAT Simulation
 Compiled 2025:10:06 with intel Linux 
 Current random number generator: ranmar
 Energy resolution as function of 1/sqrt(E)
 Header file: spect_simulation_kidney_0.h00
 Linear angle sampling within acceptance angle
 Inifile: simind.ini
 Command: spect_simulation spect_simulation_kidney_0 /fd:spect_preprocessing_atn_av.bin/fs:pbpk_kidney_act_av.bin/in:x22,3x/nn:1.027678370475769/cc:si-me/fi:lu177/02:10.909079706668853/05:10.909079706668853/08:21.82/10:53.30/14:-7/15:-7/28:0.5/29:64/31:0.23976001739501954/34:91/42:-10/76:128/77:43.636318826675414/78:256/79:256/rr:0
