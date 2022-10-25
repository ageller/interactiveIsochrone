# Interactive Isochrone
Testing an interactive plotting tool where the user can move an isochrone around to match by eye photometry data.  This is intended to adjust the starting point, and possibly also the variance on the priors, in a more sophisticated fitting routine, e.g., using BASE-9. 

Before running this, you need to download a model.  Available BASE-9 models are located here: [https://github.com/BayesianStellarEvolution/base-models](https://github.com/BayesianStellarEvolution/base-models).

To run this in a Jupyter notebook (intended purpose):

```
layout = createInteractiveIsochrone("filename.phot","isochrone.model")

def bkapp(doc):
	doc.add_root(layout)

show(bkapp)
```
