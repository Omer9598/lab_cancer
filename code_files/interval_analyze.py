import os

import matplotlib.pyplot as plt


def create_intervals(haplotype_dict):
    """
    This function will create a list for each child in the following format:
    [(interval number is the index) {start position: , end position: ,
    haplotype: (1 or 2)}
    each interval starts with the position of a variant from haplotype 1 or 2,
    and ends when the next variant is from the opposite haplotype, where a new
    interval will start
    """
    intervals = []
    current_interval = None
    for position, value in haplotype_dict.items():
        chromosome = haplotype_dict[position][-3]
        cur_haplotype = value[-2]
        if current_interval is None:
            # Start a new interval
            current_interval = {"start": position, "end": position,
                                "haplotype": cur_haplotype}
        elif cur_haplotype == current_interval["haplotype"]:
            # Continue the current interval
            current_interval["end"] = position
        else:
            # Start a new interval as haplotype changed
            intervals.append({"start": current_interval["start"],
                              "end": current_interval["end"],
                              "haplotype": current_interval["haplotype"],
                              "chromosome": chromosome})
            current_interval = {"start": position, "end": position,
                                "haplotype": cur_haplotype}
    # Add the last interval
    if current_interval is not None:
        intervals.append({"start": current_interval["start"],
                          "end": current_interval["end"],
                          "haplotype": current_interval["haplotype"],
                          "chromosome": chromosome})
    return intervals


def shared_interval(interval_lists):
    """
    This function creates a new list containing intervals that are shared
    in all the lists given, according to the haplotype.
    It also adds certainty level of 1 to the shared intervals.
    If haplotypes are different, haplotype is set to 0 and certainty level to -1.
    """
    # Initialize shared_intervals with the intervals from the first list
    shared_intervals = interval_lists[0]

    # Iterate through the remaining lists
    for interval_list in interval_lists[1:]:
        # Use a list comprehension to find shared intervals
        shared_intervals = [
            {
                "start": max(interval_1["start"], interval_2["start"]),
                "end": min(interval_1["end"], interval_2["end"]),
                "haplotype": 0 if interval_1["haplotype"] != interval_2["haplotype"] else interval_1["haplotype"],
                "chromosome": interval_1["chromosome"],
                "certainty_level": -1 if interval_1["haplotype"] != interval_2["haplotype"] else 1
            }
            for interval_1 in shared_intervals
            for interval_2 in interval_list
            if interval_1["start"] <= interval_2["end"] and interval_1["end"] >= interval_2["start"]
        ]

    return shared_intervals


def plot_interval(interval_list, plot_title, save_dir):
    """
    This function plots intervals as straight lines using Matplotlib.
    """
    fig, ax = plt.subplots()
    for i, interval in enumerate(interval_list):
        start_position = interval["start"]
        end_position = interval["end"]
        haplotype = interval["haplotype"]
        certainty_level = interval.get("certainty_level", None)
        # Determine color based on certainty level
        color = 'red' if certainty_level == -1 else 'green'
        # Plot lines for each interval
        ax.plot([start_position, end_position], [haplotype, haplotype], label=f'Interval {i + 1}', color=color)
    ax.set(xlabel='Chromosome Position', ylabel='Haplotype', title=plot_title)
    # Save the plot
    save_path = os.path.join(save_dir, f'{plot_title.replace(" ", "_")}_plot.png')
    # Create the directory if it doesn't exist
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(save_path)
