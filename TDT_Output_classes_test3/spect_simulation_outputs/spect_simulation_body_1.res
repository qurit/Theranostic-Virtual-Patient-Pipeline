


              SIMIND Monte Carlo Simulation Program    V8.0  
------------------------------------------------------------------------------
 Phantom S : h2o       Crystal...: nai       InputFile.: spect_simulation  
 Phantom B : h2o       BackScatt.: pmt       OutputFile: spect_simulation_b
 Collimator: pb_sb2    SourceRout: smap      SourceImg.: pbpk_body_act_av  
 Cover.....: al        ScoreRout.: scattwin  DensityImg: spect_preprocessin
------------------------------------------------------------------------------
 PhotonEnergy.......: 200          lu177     PhotonsPerProj....: 178317         
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
 Emission type......: 2                      Initial Weight....: 1.26964        
 NN ScalingFactor...: 8.71286                Energy Channels...: 512            
                                                                              
 SPECT DATA
 RotationMode.......: -360                   Nr of Projections.: 64             
 RotationAngle......: 5.625                  Projection.[start]: 1              
 Orbital fraction...: 0                      Projection...[end]: 64             
 Center of Rotation File: spect_simulation_body_1.cor
                                                                              
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
   1   0.404E+02 0.356E+02 0.473E+01 0.754E+01 0.883E+00 0.631E+00
   2   0.233E+03 0.439E+02 0.189E+03 0.232E+00 0.188E+00 0.364E+01
   3   0.283E+01 0.766E+00 0.206E+01 0.372E+00 0.271E+00 0.442E-01
  
  Win  Geo(Air)  Pen(Air)  Sca(Air)  Geo(Tot)  Pen(Tot)  Sca(Tot)
   1   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
   2   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
   3   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
  
  Win   SC 1  SC 2  SC 3
   1   72.0% 23.8%  4.1%
   2   88.3% 10.7%  1.0%
   3   79.3% 17.3%  3.5%
                                                                              
 INTERACTIONS IN THE CRYSTAL
 MaxValue spectrum..: 10.73          
 MaxValue projection: 0.5526E-02     
 CountRate spectrum.: 19.86          
 CountRate E-Window.: 3.789          
                                                                              
 SCATTER IN ENERGY WINDOW
 Scatter/Primary....: 0.31173        
 Scatter/Total......: 0.23765        
 Scatter order 1....: 85.38 %        
 Scatter order 2....: 13.18 %        
 Scatter order 3....: 1.44 %         
                                                                              
 CALCULATED DETECTOR PARAMETERS
 Efficiency E-window: 0.1705         
 Efficiency spectrum: 0.8938         
 Sensitivity Cps/MBq: 3.7887         
 Sensitivity Cpm/uCi: 8.4109         
                                                                              
 Simulation started.: 2025:12:13 00:33:44
 Simulation stopped.: 2025:12:13 00:36:13
 Elapsed time.......: 0 h, 2 m and 29 s
 DetectorHits.......: 385437         
 DetectorHits/CPUsec: 3500           
                                                                              
 OTHER INFORMATION
 XCAT Simulation
 Compiled 2025:10:06 with intel Linux 
 Current random number generator: ranmar
 Energy resolution as function of 1/sqrt(E)
 Header file: spect_simulation_body_1.h00
 Linear angle sampling within acceptance angle
 Inifile: simind.ini
 Command: spect_simulation spect_simulation_body_1 /fd:spect_preprocessing_atn_av.bin/fs:pbpk_body_act_av.bin/in:x22,3x/nn:8.71286392211914/cc:si-me/fi:lu177/02:10.909079706668853/05:10.909079706668853/08:21.82/10:53.30/14:-7/15:-7/28:0.5/29:64/31:0.23976001739501954/34:91/42:-10/76:128/77:43.636318826675414/78:256/79:256/rr:1
