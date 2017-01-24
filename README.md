# mangum_galaxies

Script to derive Moment0, Moment1, and Moment2 from a set of
input-defined spectral lines in an image cube.  Currently simply
calculates moments over a defined HWZI for each line in band. 

To run in ipython use:

run ~/Python/CubeLineMoment.py




notes:

    # From NGC4945 H2COJ32K02 spectral baseline
    #noisemap = cube[1:50,:,:].std(axis=0)
    # For NGC253 H2COJ54K23...
    #noisemap = cube[35:45,:,:].std(axis=0)
    # For NGC253 H2COJ54K1...
    #noisemap = cube[50:90,:,:].std(axis=0)
    #
    # Note that "line widths" in the following are HWZI values...
    #
    # NGC253 H2CO J=3-2 K=0,2
    #
    #    my_line_list = [217.289800, 217.299162, 217.467150, 217.517110, 217.802057, 217.88639, 217.943821, 218.15897, 218.222192, 218.324711, 218.440050, 218.475632, 218.760071, 218.85439, 218.9033555, 218.981019] * u.GHz
    #    my_line_widths = [50.0, 50.0, 60.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 50.0, 40.0, 40.0] * u.km/u.s
    #    my_line_names = ['13CNF122','CH3OH67','13CNF132','CH3OCHO88','CH3OCHO4847','CH3OH2020','CH3OCHO4546','CH3OCHO??','H2COJ32K0','HC3N2423v0','CH3OH43','H2COJ32K221','H2COJ32K210','HC3N2423v6','OCS1817','HNCO109']
    # These are:
    # 13CN N=2-1,J=5/2-3/2,F1=2-2,F=3-2 at 217289.800 MHz              (Blend with CH3OH 6(15)-7(26))
    # CH3OH 6(15)-7(26)            217299.162 MHz (+1268.05 km/s)      (Blend with 13CN F1=2-2)
    # 13CN N=2-1,J=5/2-3/2,F1=3-2,F=3-2  at 217467.150 MHz             (Set width larger to encompass multiplet)
    # CH3OCHO 8(5,4)-8(3,5)A (Methyl Formate) at 217517.110 MHz
    # CH3OCHO 48(14,35)-47(15,32)A (Methyl Formate) at 217802.057 MHz
    # CH3OH 20(1,19)-20(0,20) at 217886.39 MHz
    # CH3OCHO 45(29,16)-46(28,18)E (Methyl Formate) 217943.821 MHz
    # CH3OCHO ?? (Methyl Formate) at 218158.97 MHz
    # H2CO 3(03)-2(02)             218222.192 MHz (0.0 km/s)           yes
    # HC3N 24-23                   218324.711 MHz (-140.84 km/s)       yes
    # CH3OH 4(22)-3(12)E           218440.050 MHz (-304.79 km/s)       yes (blend with H2CO K2)
    # H2CO 3(22)-2(21)             218475.632 MHz (-348.17 km/s)       yes
    # H2CO 3(21)-2(20)             218760.071 MHz (-738.94 km/s)       yes
    # HC3N 24-23 v6=1              218854.39 MHz                       (Blend with v7=1; stronger of two vib lines?)
    # HC3N 24-23 v7=1              218860.80 MHz                       (Blend with v6=1)
    # OCS 18-17 at 218903.3555 MHz
    # HNCO 10(1,10)-9(1,9) at 218981.019 MHz
    #
    #
    # NGC4945 H2CO J=3-2 K=0,2
    #
    #my_line_list = [217.299162, 217.802057, 217.467150, 217.296605, 217.943821, 218.222192, 218.324711, 218.440050, 218.475632, 218.760071] * u.GHz
    #my_line_widths = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0] * u.km/u.s
    #my_line_names = ['CH3OH67','CH3OCHO4847','13CN21-1','13CN21-2','CH3OCHO4546','H2COJ32K0','HC3N2423','CH3OH43','H2COJ32K221','H2COJ32K210']
    # These are:
    # 13CN N=2-1 (2)               217296.605 MHz (+1732 km/s)        yes
    # CH3OH 6(15)-7(26)            217299.162 MHz (+1268.05 km/s)     yes
    # 13CN N=2-1 (1)               217467.150 MHz (+1500 km/s)        yes
    # CH3OCHO 48(14,35)-47(15,32)A 217802.057 MHz (+1045 km/s)        yes
    # CH3OCHO 45(29,16)-46(28,18)E 217943.821 MHz (+382 km/s)         yes
    # H2CO 3(03)-2(02)             218222.192 MHz (0.0 km/s)          yes
    # HC3N 24-23                   218324.711 MHz (-140.84 km/s)      yes
    # CH3OH 4(22)-3(12)E           218440.050 MHz (-304.79 km/s)      yes (blend with H2CO K2)
    # H2CO 3(22)-2(21)             218475.632 MHz (-348.17 km/s)      yes
    # H2CO 3(21)-2(20)             218760.071 MHz (-738.94 km/s)      yes
    #
    #
    # H213CO J=3-2 K=1
    #
    #my_line_list = [220.398684, 219.908486, 219.798274, 219.675114, 219.560358] * u.GHz
    #my_line_widths = [150.0, 80.0, 80.0, 80.0, 80.0] * u.km/u.s
    #my_line_names = ['13CO21','H213COJ32K1','HNCO109','HC3N2423','C18O21']
    # These are:
    # 13CO 2-1                       220398.684 MHz                     yes
    # H213CO 3(12)-2(11)             219908.486 MHz                     yes
    # HNCO 10(0,10)-9(0,9)           219798.274 MHz                     yes
    # HC3N 24-23                     219675.114 MHz                     yes
    # C180 2-1                       219560.358 MHz                     yes
    #
    #
    # H2CO J=5-4 K=1
    #
    #my_line_list = [352.199050, 351.984515, 352.603542, 351.768645, 351.6334, 351.45424, 351.236343, 351.047000, 350.905070] * u.GHz
    #my_line_widths = [40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0] * u.km/u.s
    #my_line_names = ['CH3OCHO4545','CH3OCHO3231','CH3OCHO3029','H2COJ54K1','HNCO1615','CH2NH1010','CH3OH910','Unidentified351236','CH3OH10']
    # These are:
    # CH3OCHO 45(14,32)-45(13,33)E 352199.050 MHz                     yes
    # CH3OCHO 32(1,31)-31(1,30)E   351984.515 MHz                     yes
    # CH3OCHO 30(4,27)-29(3,26)A   352603.542 MHz                     yes
    # H2CO 5(15)-4(14)             351768.645 MHz                     yes
    # HNCO 16(0,16)-15(0,15)       351633.4 MHz (+115.264 km/s)       yes (next to H2CO)              2e14
    # CH2NH 10(1,9)-10(0,10)       351454.24 MHz (+267.95 km/s)       yes (some positions)
    # CH3OH 9(5,5)-10(4,6) E       351236.343 MHz (+453.65 km/s)      yes (some positions)            1e15
    # UNIDENTIFIED                 351047.000 MHz (+615.02 km/s)      yes (strong line next to
    #                                                                      CH3OH 1-0)
    # CH3OH 1(1,1)-0(0,0) A++      350905.070 MHz (+735.98 km/s)      yes
