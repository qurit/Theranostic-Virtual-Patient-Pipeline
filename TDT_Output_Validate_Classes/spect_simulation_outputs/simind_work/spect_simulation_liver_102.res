


              SIMIND Monte Carlo Simulation Program    V8.0  
------------------------------------------------------------------------------
 Phantom S : h2o       Crystal...: nai       InputFile.: spect_simulation  
 Phantom B : h2o       BackScatt.: pmt       OutputFile: spect_simulation_l
 Collimator: pb_sb2    SourceRout: smap      SourceImg.: pbpk_liver_act_av 
 Cover.....: al        ScoreRout.: scattwin  DensityImg: spect_preprocessin
------------------------------------------------------------------------------
 PhotonEnergy.......: 200          lu177     PhotonsPerProj....: 531295         
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
 Emission type......: 2                      Initial Weight....: 0.42612        
 NN ScalingFactor...: 91.25652               Energy Channels...: 512            
                                                                              
 SPECT DATA
 RotationMode.......: -360                   Nr of Projections.: 64             
 RotationAngle......: 5.625                  Projection.[start]: 1              
 Orbital fraction...: 0                      Projection...[end]: 64             
 Center of Rotation File: spect_simulation_liver_102.cor
                                                                              
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
   1   0.553E+02 0.526E+02 0.274E+01 0.192E+02 0.950E+00 0.864E+00
   2   0.165E+03 0.565E+02 0.109E+03 0.521E+00 0.342E+00 0.258E+01
   3   0.235E+01 0.109E+01 0.125E+01 0.875E+00 0.467E+00 0.366E-01
  
  Win  Geo(Air)  Pen(Air)  Sca(Air)  Geo(Tot)  Pen(Tot)  Sca(Tot)
   1   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
   2   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
   3   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
  
  Win   SC 1  SC 2  SC 3
   1   63.6% 30.0%  6.4%
   2   84.1% 14.2%  1.7%
   3   72.5% 22.7%  4.8%
                                                                              
 INTERACTIONS IN THE CRYSTAL
 MaxValue spectrum..: 7.748          
 MaxValue projection: 0.7330E-02     
 CountRate spectrum.: 18.54          
 CountRate E-Window.: 2.847          
                                                                              
 SCATTER IN ENERGY WINDOW
 Scatter/Primary....: 0.71944        
 Scatter/Total......: 0.41842        
 Scatter order 1....: 80.18 %        
 Scatter order 2....: 17.38 %        
 Scatter order 3....: 2.44 %         
                                                                              
 CALCULATED DETECTOR PARAMETERS
 Efficiency E-window: 0.1394         
 Efficiency spectrum: 0.908          
 Sensitivity Cps/MBq: 2.8467         
 Sensitivity Cpm/uCi: 6.3197         
                                                                              
 Simulation started.: 2025:12:19 15:17:53
 Simulation stopped.: 2025:12:19 15:23:58
 Elapsed time.......: 0 h, 6 m and 5 s
 DetectorHits.......: 1721338        
 DetectorHits/CPUsec: 4745           
                                                                              
 OTHER INFORMATION
 XCAT Simulation
 Compiled 2025:10:06 with intel Linux 
 Current random number generator: ranmar
 Energy resolution as function of 1/sqrt(E)
 Header file: spect_simulation_liver_102.h00
 Linear angle sampling within acceptance angle
 Inifile: simind.ini
 Command: spect_simulation spect_simulation_liver_102 /fd:spect_preprocessing_atn_av.bin/fs:pbpk_liver_act_av.bin/in:x22,3x/nn:91.25652313232422/cc:si-me/fi:lu177/02:55.7617868158035/05:55.7617868158035/08:111.52/10:53.30/14:-7/15:-7/28:0.5/29:64/31:0.239760018279776/34:91/42:-10/76:128/77:223.047147263214/78:256/79:256/rr:102
