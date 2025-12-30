


              SIMIND Monte Carlo Simulation Program    V8.0  
------------------------------------------------------------------------------
 Phantom S : h2o       Crystal...: nai       InputFile.: spect_simulation  
 Phantom B : h2o       BackScatt.: pmt       OutputFile: spect_simulation_b
 Collimator: pb_sb2    SourceRout: smap      SourceImg.: pbpk_body_act_av  
 Cover.....: al        ScoreRout.: scattwin  DensityImg: spect_preprocessin
------------------------------------------------------------------------------
 PhotonEnergy.......: 200          lu177     PhotonsPerProj....: 6368480        
 EnergyResolution...: 8.8          Spectra   Activity..........: 1              
 MaxScatterOrder....: 3            si-me     DetectorLenght....: 111.52         
 DetectorWidth......: 53.3         SPECT     DetectorHeight....: 0.953          
 UpperEneWindowTresh: 220          x-rays    Distance to det...: 35.203         
 LowerEneWindowTresh: 180          BScatt    ShiftSource X.....: 0              
 PixelSize  I.......: 0.5          Random    ShiftSource Y.....: 0              
 PixelSize  J.......: 0.5          Cover     ShiftSource Z.....: 0              
 HalfLength S.......: 55.76179     Phantom   HalfLength P......: 55.76179       
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
 MatrixSize J.......: 223                    AcceptanceAngle...: 4.13771        
 Emission type......: 2                      Initial Weight....: 0.03555        
 NN ScalingFactor...: 311.17368              Energy Channels...: 512            
                                                                              
 SPECT DATA
 RotationMode.......: -360                   Nr of Projections.: 64             
 RotationAngle......: 5.625                  Projection.[start]: 1              
 Orbital fraction...: 0                      Projection...[end]: 64             
 Center of Rotation File: spect_simulation_body_21.cor
                                                                              
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
 CT-Pixel size......: 0.23976                Slice thickness...: 1.22553        
 StartImage.........: 1                      No of CT-Images...: 91             
 MatrixSize I.......: 256                    CTmapOrientation..: 0              
 MatrixSize J.......: 256                    StepSize..........: 0.1            
 CenterPoint I......: 129                    ShiftPhantom X....: 0              
 CenterPoint J......: 129                    ShiftPhantom Y....: 0              
 CenterPoint K......: 46.5                   ShiftPhantom Z....: 0              
                                                                              
 INFO FOR TCT file
 MatrixSize I.......: 128                    MatrixSize J......: 128            
 MatrixSize K.......: 223                    Units.............: mu                  
 Scout File.........: F
                                                                              
------------------------------------------------------------------------------
  Scattwin results: Window file: spect_simulation.win
  
  Win  WinAdded  Range(keV)   ScaleFactor
   1       0    169.4 - 187.2   1.000
   2       0    187.2 - 228.8   1.000
   3       0    228.8 - 252.9   1.000
  
  Win    Total    Scatter   Primary  S/P-Ratio S/T Ratio  Cps/MBq
   1   0.549E+02 0.497E+02 0.521E+01 0.953E+01 0.905E+00 0.857E+00
   2   0.262E+03 0.542E+02 0.208E+03 0.260E+00 0.207E+00 0.410E+01
   3   0.332E+01 0.101E+01 0.230E+01 0.440E+00 0.305E+00 0.518E-01
  
  Win  Geo(Air)  Pen(Air)  Sca(Air)  Geo(Tot)  Pen(Tot)  Sca(Tot)
   1   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
   2   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
   3   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
  
  Win   SC 1  SC 2  SC 3
   1   69.8% 25.3%  4.9%
   2   87.1% 11.6%  1.3%
   3   77.0% 19.1%  3.9%
                                                                              
 INTERACTIONS IN THE CRYSTAL
 MaxValue spectrum..: 12.21          
 MaxValue projection: 0.7477E-03     
 CountRate spectrum.: 24.72          
 CountRate E-Window.: 4.312          
                                                                              
 SCATTER IN ENERGY WINDOW
 Scatter/Primary....: 0.35841        
 Scatter/Total......: 0.26385        
 Scatter order 1....: 83.84 %        
 Scatter order 2....: 14.31 %        
 Scatter order 3....: 1.85 %         
                                                                              
 CALCULATED DETECTOR PARAMETERS
 Efficiency E-window: 0.1566         
 Efficiency spectrum: 0.8978         
 Sensitivity Cps/MBq: 4.312          
 Sensitivity Cpm/uCi: 9.5726         
                                                                              
 Simulation started.: 2025:12:19 15:25:06
 Simulation stopped.: 2025:12:19 16:21:12
 Elapsed time.......: 0 h, 56 m and 6 s
 DetectorHits.......: 17519743       
 DetectorHits/CPUsec: 5232           
                                                                              
 OTHER INFORMATION
 XCAT Simulation
 Compiled 2025:10:06 with intel Linux 
 Current random number generator: ranmar
 Energy resolution as function of 1/sqrt(E)
 Header file: spect_simulation_body_21.h00
 Linear angle sampling within acceptance angle
 Inifile: simind.ini
 Command: spect_simulation spect_simulation_body_21 /fd:spect_preprocessing_atn_av.bin/fs:pbpk_body_act_av.bin/in:x22,3x/nn:311.1736755371094/cc:si-me/fi:lu177/02:55.7617868158035/05:55.7617868158035/08:111.52/10:53.30/14:-7/15:-7/28:0.5/29:64/31:0.239760018279776/34:91/42:-10/76:128/77:223.047147263214/78:256/79:256/rr:21
