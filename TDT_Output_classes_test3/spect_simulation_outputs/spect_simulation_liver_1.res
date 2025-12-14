


              SIMIND Monte Carlo Simulation Program    V8.0  
------------------------------------------------------------------------------
 Phantom S : h2o       Crystal...: nai       InputFile.: spect_simulation  
 Phantom B : h2o       BackScatt.: pmt       OutputFile: spect_simulation_l
 Collimator: pb_sb2    SourceRout: smap      SourceImg.: pbpk_liver_act_av 
 Cover.....: al        ScoreRout.: scattwin  DensityImg: spect_preprocessin
------------------------------------------------------------------------------
 PhotonEnergy.......: 200          lu177     PhotonsPerProj....: 14876          
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
 Emission type......: 2                      Initial Weight....: 15.21901       
 NN ScalingFactor...: 2.55518                Energy Channels...: 512            
                                                                              
 SPECT DATA
 RotationMode.......: -360                   Nr of Projections.: 64             
 RotationAngle......: 5.625                  Projection.[start]: 1              
 Orbital fraction...: 0                      Projection...[end]: 64             
 Center of Rotation File: spect_simulation_liver_1.cor
                                                                              
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
   1   0.368E+02 0.346E+02 0.226E+01 0.153E+02 0.939E+00 0.576E+00
   2   0.130E+03 0.419E+02 0.883E+02 0.475E+00 0.322E+00 0.204E+01
   3   0.187E+01 0.786E+00 0.108E+01 0.729E+00 0.421E+00 0.292E-01
  
  Win  Geo(Air)  Pen(Air)  Sca(Air)  Geo(Tot)  Pen(Tot)  Sca(Tot)
   1   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
   2   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
   3   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
  
  Win   SC 1  SC 2  SC 3
   1   68.3% 27.1%  4.6%
   2   86.3% 12.5%  1.1%
   3   75.1% 23.5%  1.4%
                                                                              
 INTERACTIONS IN THE CRYSTAL
 MaxValue spectrum..: 5.815          
 MaxValue projection: 0.2481E-01     
 CountRate spectrum.: 13.32          
 CountRate E-Window.: 2.210          
                                                                              
 SCATTER IN ENERGY WINDOW
 Scatter/Primary....: 0.63869        
 Scatter/Total......: 0.38976        
 Scatter order 1....: 83.09 %        
 Scatter order 2....: 15.23 %        
 Scatter order 3....: 1.67 %         
                                                                              
 CALCULATED DETECTOR PARAMETERS
 Efficiency E-window: 0.1501         
 Efficiency spectrum: 0.9049         
 Sensitivity Cps/MBq: 2.2097         
 Sensitivity Cpm/uCi: 4.9055         
                                                                              
 Simulation started.: 2025:12:13 00:40:50
 Simulation stopped.: 2025:12:13 00:41:08
 Elapsed time.......: 0 h, 0 m and 18 s
 DetectorHits.......: 32240          
 DetectorHits/CPUsec: 2184           
                                                                              
 OTHER INFORMATION
 XCAT Simulation
 Compiled 2025:10:06 with intel Linux 
 Current random number generator: ranmar
 Energy resolution as function of 1/sqrt(E)
 Header file: spect_simulation_liver_1.h00
 Linear angle sampling within acceptance angle
 Inifile: simind.ini
 Command: spect_simulation spect_simulation_liver_1 /fd:spect_preprocessing_atn_av.bin/fs:pbpk_liver_act_av.bin/in:x22,3x/nn:2.555182695388794/cc:si-me/fi:lu177/02:10.909079706668853/05:10.909079706668853/08:21.82/10:53.30/14:-7/15:-7/28:0.5/29:64/31:0.23976001739501954/34:91/42:-10/76:128/77:43.636318826675414/78:256/79:256/rr:1
