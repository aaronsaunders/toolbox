def plot_length(fname, sub=False):
    """
    Parses a sequence file and returns a plot of sequence lengths.
    Optional argument to subset the file.
    
    """
    import matplotlib.pyplot as plt
    from toolbox.bioutils import get_sizes
    
    sizes = get_sizes(fname, sub)

    plt.hist(sizes)
    plt.title("%s\n%i sequences, range %i to %i bp" \
                % (fname, len(sizes), min(sizes),max(sizes)))
    plt.xlabel("Sequence length (bp)")
    plt.ylabel("Count")
    plt.show()

    return
