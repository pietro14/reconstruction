Giorgio was reconstructing iron55 with chanvese:
Paramenters:
	Configfile:
		line 36: 'nsigma' = 2
		line 53:  'vector_eps' = [1, 6.5,  4,  2]
	supercluster.py:
		line 52: gimage =inverse_gaussian_gradient(clustered_data,5,3)
		line 66: ls = morphological_chan_vese(gimage, 400, init_ls,
                                     smoothing=1, lambda1=2,
                                     lambda2=1.5)

Samuele was reconstructing for directionality:
Paramenters:
	Configfile:
		line 36: 'nsigma' = 1.5
		line 53:  'vector_eps' = [1, 3.5,  6.5,  3.5]
	supercluster.py:
		line 52: gimage =inverse_gaussian_gradient(clustered_data,8,1.4)
		line 66: ls = morphological_chan_vese(gimage, 400, init_ls,
                                     smoothing=1, lambda1=2,
                                     lambda2=1.5)

Flaminia was reconstructing nuclear recoils:
Paramenters:
	Configfile:
		line 36: 'nsigma' = 1.5
		line 53:  'vector_eps' = [1, 3.5,  6.5,  3.5]
	supercluster.py:
		line 52: gimage =inverse_gaussian_gradient(clustered_data,8,1.4)
		line 66: ls = morphological_chan_vese(gimage, 400, init_ls,
                                     smoothing=1, lambda1=2,
                                     lambda2=1.5)

Atul was reconstructing everything:
Paramenters:
	Configfile:
		line 36: 'nsigma' = 
		line 53:  'vector_eps' = [, ,  ,  ]
	supercluster.py:
		line 52: gimage =inverse_gaussian_gradient(clustered_data,,)
		line 66: ls = morphological_chan_vese(gimage, 400, init_ls,
                                     smoothing=1, lambda1=2,
                                     lambda2=1.5)
