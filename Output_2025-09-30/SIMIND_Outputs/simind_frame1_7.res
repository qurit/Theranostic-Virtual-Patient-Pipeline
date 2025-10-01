


              SIMIND Monte Carlo Simulation Program    V8.0  
------------------------------------------------------------------------------
 Phantom S : h2o       Crystal...: nai       InputFile.: SIMIND            
 Phantom B : h2o       BackScatt.: pmt       OutputFile: simind_frame1_7   
 Collimator: pb_sb2    SourceRout: smap      SourceImg.: PBPK_200_act_av   
 Cover.....: al        ScoreRout.: scattwin  DensityImg: for_simind_atn_av 
------------------------------------------------------------------------------
 PhotonEnergy.......: 200          lu177     PhotonsPerProj....: 124999         
 EnergyResolution...: 8.8          Spectra   Activity..........: 1              
 MaxScatterOrder....: 3            si-me     DetectorLenght....: 25             
 DetectorWidth......: 0            SPECT     DetectorHeight....: 0.953          
 UpperEneWindowTresh: 220          x-rays    Distance to det...: 10             
 LowerEneWindowTresh: 180          BScatt    ShiftSource X.....: 0              
 PixelSize  I.......: 0.5          Random    ShiftSource Y.....: 0              
 PixelSize  J.......: 0.5          Cover     ShiftSource Z.....: 0              
 HalfLength S.......: 83.712       Phantom   HalfLength P......: 83.712         
 HalfWidth  S.......: 16           Resolut   HalfWidth  P......: 16             
 HalfHeight S.......: 16           Header    HalfHeight P......: 16             
 SourceType.........: XcatBinMap   SaveMap   PhantomType.......: XcatBinMap   
------------------------------------------------------------------------------
 GENERAL DATA
 keV/channel........: 0.5                    CutoffEnergy......: 0              
 Photons/Bq.........: 0.2264                 StartingAngle.....: 0              
 CameraOffset X.....: 0                      CoverThickness....: 0.1            
 CameraOffset Y.....: 0                      BackscatterThickn.: 10             
 MatrixSize I.......: 128                    IntrinsicResolut..: 0.34           
 MatrixSize J.......: 334                    AcceptanceAngle...: 4.13771        
 Emission type......: 2                      Initial Weight....: 1.8112         
 NN ScalingFactor...: 0.00284                Energy Channels...: 512            
                                                                              
 SPECT DATA
 RotationMode.......: -360                   Nr of Projections.: 10             
 RotationAngle......: 36                     Projection.[start]: 1              
 Orbital fraction...: 1                      Projection...[end]: 10             
                                                                              
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
 RotationCentre.....:  65, 65                Bone definition...: 1170           
 CT-Pixel size......: 0.39062                Slice thickness...: 1.308          
 StartImage.........: 1                      No of CT-Images...: 128            
 MatrixSize I.......: 128                    CTmapOrientation..: 0              
 MatrixSize J.......: 128                    StepSize..........: 0.1            
 CenterPoint I......: 65                     ShiftPhantom X....: 0              
 CenterPoint J......: 65                     ShiftPhantom Y....: 0              
 CenterPoint K......: 65                     ShiftPhantom Z....: 0              
                                                                              
 INFO FOR TCT file
 MatrixSize I.......: 128                    MatrixSize J......: 128            
 MatrixSize K.......: 334                    Units.............: HU                  
 Scout File.........: F
 CT unit number.....: 4                      CT Voltage......kV: 120            
 CT effective hv....: 58.1           
                                                                              
------------------------------------------------------------------------------
  Scattwin results: Window file: SIMIND.win          
  
  Win  WinAdded  Range(keV)   ScaleFactor
   1       0    169.4 - 187.2   1.000
   2       0    187.2 - 228.8   1.000
   3       0    228.8 - 252.9   1.000
  
  Win    Total    Scatter   Primary  S/P-Ratio S/T Ratio  Cps/MBq
   1   0.000E+00 0.000E+00 0.000E+00       NaN       NaN       NaN
   2   0.000E+00 0.000E+00 0.000E+00       NaN       NaN       NaN
   3   0.000E+00 0.000E+00 0.000E+00       NaN       NaN       NaN
  
  Win  Geo(Air)  Pen(Air)  Sca(Air)  Geo(Tot)  Pen(Tot)  Sca(Tot)
   1     0.00%     0.00%     0.00%     0.00%     0.00%     0.00%
   2     0.00%     0.00%     0.00%     0.00%     0.00%     0.00%
   3     0.00%     0.00%     0.00%     0.00%     0.00%     0.00%
  
                                                                              
 INTERACTIONS IN THE CRYSTAL
 MaxValue spectrum..: 0.000          
 MaxValue projection: 0.000          
 CountRate spectrum.: 0.000          
 CountRate E-Window.: 0.000          
                                                                              
 CALCULATED DETECTOR PARAMETERS
 Efficiency E-window: 0              
 Efficiency spectrum: 0              
 Sensitivity Cps/MBq: 0              
 Sensitivity Cpm/uCi: 0              
                                                                              
 Simulation started.: 2025:09:30 22:41:24
 Simulation stopped.: 2025:09:30 22:41:32
 Elapsed time.......: 0 h, 0 m and 8 s
 DetectorHits.......: 0              
 DetectorHits/CPUsec: 0              
                                                                              
 OTHER INFORMATION
 EMISSION
 Compiled 2025:01:28 with INTEL Mac   
 Current random number generator: ranmar
 Energy resolution as function of 1/sqrt(E)
 Header file: simind_frame1_7.h00
 Linear angle sampling within acceptance angle
 Inifile: simind.ini
 Command: SIMIND SIMIND_frame1_7 /fd:for_simind_atn_av.bin/fs:PBPK_200_act_av.bin/14:-7/15:-7/cc:si-me/fi:lu177/nn:0.0028433240950107574/in:x22,6x/02:83.71199951171876/05:83.71199951171876/12:10/28:0.5/29:10/31:0.39062480926513676/32:0/34:128/76:128/77:334.84799804687503/78:128/79:128
