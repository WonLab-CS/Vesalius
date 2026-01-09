# The vesalius_assay class

The vesalius_assay class is the functional unit of vesalius. Each assay
is stored within this class and it contains all the required information
to run analysis on your assay of choice. In this object, you can find
spatial tiles, image embeddings, spatial territories, differentially
expressed genes (DEG), count matrices (raw and normalised), microscopy
images (if present) and a functional log that lets you see what had been
run on this object.

## Slots

- `assay`:

  character assay name

- `tiles`:

  data.frame containing spatial coordinates and pixels tiles once they
  have been computed

- `embeddings`:

  list containing latent space embeddings in the form of data.frames.

- `active`:

  matrix containing active embedding data

- `territories`:

  data.frame containing spatial color segments, spatial territories, or
  layers.

- `DEG`:

  list of data.frame for each differentially gene expression trial

- `counts`:

  list that containing count matrices. Raw and normalised will be stored
  here and named by the normalisation method used.

- `meta`:

  list containing metadata associated with the assay

- `cost`:

  list containing cost matrices used for mapping assays

- `map`:

  data.frame containing mapping results between assays

- `image`:

  list containing associated microscopy images (NOT implemented)

- `log`:

  list containing analysis history of the object.
