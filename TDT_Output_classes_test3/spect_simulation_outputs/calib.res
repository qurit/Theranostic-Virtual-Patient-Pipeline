


              SIMIND Monte Carlo Simulation Program    V8.0  
------------------------------------------------------------------------------
 Phantom S : h2o       Crystal...: nai       InputFile.: jaszak            
 Phantom B : h2o       BackScatt.: lucite    OutputFile: calib             
 Collimator: pb_sb2    SourceRout: none      SourceImg.: none              
 Cover.....: al        ScoreRout.: scattwin  DensityImg: none              
------------------------------------------------------------------------------
 PhotonEnergy.......: 208          lu177     PhotonsPerProj....: 20000000       
 EnergyResolution...: 10           si-me     Activity..........: 1              
 MaxScatterOrder....: 0            SPECT     DetectorLenght....: 25             
 DetectorWidth......: 0            Random    DetectorHeight....: 0.9525         
 UpperEneWindowTresh: 228.8        Resolut   Distance to det...: 15             
 LowerEneWindowTresh: 187.2                  ShiftSource X.....: 0              
 PixelSize  I.......: 0.4                    ShiftSource Y.....: 0              
 PixelSize  J.......: 0.4                    ShiftSource Z.....: 0              
 HalfLength S.......: 0                      HalfLength P......: 10             
 HalfWidth  S.......: 0                      HalfWidth  P......: 11             
 HalfHeight S.......: 0                      HalfHeight P......: 11             
 SourceType.........: PointSrc               PhantomType.......: HorCylinder  
------------------------------------------------------------------------------
 GENERAL DATA
 keV/channel........: 1                      CutoffEnergy......: 0              
 Photons/Bq.........: 0.2264                 StartingAngle.....: 0              
 CameraOffset X.....: 0                      CoverThickness....: 0              
 CameraOffset Y.....: 0                      BackscatterThickn.: 0              
 MatrixSize I.......: 128                    IntrinsicResolut..: 0.36           
 MatrixSize J.......: 128                    AcceptanceAngle...: 4.13771        
 Emission type......: 2                      Initial Weight....: 0.01132        
 NN ScalingFactor...: 1                      Energy Channels...: 512            
                                                                              
 SPECT DATA
 RotationMode.......: -360                   Nr of Projections.: 1              
 RotationAngle......: 360                    Projection.[start]: 1              
 Orbital fraction...: 1                      Projection...[end]: 1              
                                                                              
 COLLIMATOR DATA FOR ROUTINE: Analytical          
 CollimatorCode.....: si-me                  CollimatorType....: Parallel 
 HoleSize X.........: 0.294                  Distance X........: 0.114          
 HoleSize Y.........: 0.33948                Distance Y........: 0.26847        
 CenterShift X......: 0.204                  X-Ray flag........: F              
 CenterShift Y......: 0.35334                CollimThickness...: 4.064          
 HoleShape..........: Hexagonal              Space Coll2Det....: 0              
 CollDepValue [57]..: 0                      CollDepValue [58].: 0              
 CollDepValue [59]..: 0                      CollDepValue [60].: 0              
                                                                              
------------------------------------------------------------------------------
  Scattwin results: Window file: jaszak.win          
  
  Win  WinAdded  Range(keV)   ScaleFactor
   1       0    187.2 - 228.8   1.000
  
  Win    Total    Scatter   Primary  S/P-Ratio S/T Ratio  Cps/MBq
   1   0.107E+02 0.000E+00 0.107E+02 0.000E+00 0.000E+00 0.107E+02
  
  Win  Geo(Air)  Pen(Air)  Sca(Air)  Geo(Tot)  Pen(Tot)  Sca(Tot)
   1   100.00%     0.00%     0.00%   100.00%     0.00%     0.00%
  
                                                                              
 INTERACTIONS IN THE CRYSTAL
 MaxValue spectrum..: 0.8655         
 MaxValue projection: 0.6768         
 CountRate spectrum.: 34.22          
 CountRate E-Window.: 10.66          
                                                                              
 CALCULATED DETECTOR PARAMETERS
 Efficiency E-window: 0.2525         
 Efficiency spectrum: 0.8105         
 Sensitivity Cps/MBq: 10.6609        
 Sensitivity Cpm/uCi: 23.6672        
                                                                              
 Simulation started.: 2025:12:13 00:36:14
 Simulation stopped.: 2025:12:13 00:36:25
 Elapsed time.......: 0 h, 0 m and 11 s
 DetectorHits.......: 20000000       
 DetectorHits/CPUsec: 1856931        
                                                                              
 OTHER INFORMATION
 
 Compiled 2025:10:06 with intel Linux 
 Current random number generator: ranmar
 Energy resolution as function of 1/sqrt(E)
 Linear angle sampling within acceptance angle
 Inifile: simind.ini
 Command: jaszak calib/fi:lu177/cc:si-me/29:1/15:5/fa:11/fa:15/fa:14
