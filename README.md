#A neural network model of attention

The objective is to use neural network to explain attention-induced change in neuronal RF size, position and gain.

* nnattend.m : The main function to implement simulation 
* RZutil, knkutil: two utility function repositories. Please add to your matlab path.
* extractrealvxs: RZ uses it to read, process and save vxs data. Can be ignored.
* plotcircle.m: plot a circle, useful when plotting RF.
* nnattenddata.mat: data with four field:
    * beta: 50(nNenuron) x 50(25 stim x 2 fix/att) x 100 bootstrap samples
    * betamn: take the median(beta,3), to obtain the median of responses
    * image: 800x800x25 stimuli images
    * pRFparams: a cell array with RF parameters for the fixation task in element{1}, and attention task in element{2}. Each task has a 5(params)x50(nNeurons) matrix. 5 params are (x,y),size,gain,eccentriciy. All params are in pixel unit, except gain.
