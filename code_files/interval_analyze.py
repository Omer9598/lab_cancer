import os
from collections import defaultdict

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
    This function will create a new list containing intervals that are shared
    in all the lists given, according to the haplotype.
    The function will also add certainty level of 1 to the shared intervals
    """
    # Initialize shared_intervals with the intervals from the first list
    shared_intervals = interval_lists[0]
    # Iterate through the remaining lists
    for interval_list in interval_lists[1:]:
        # A temporary list to store shared intervals for the current list
        temp_shared_intervals = []
        # Iterate through each interval in the current list
        for interval_1 in shared_intervals:
            for interval_2 in interval_list:
                if (
                        interval_1["haplotype"] == interval_2["haplotype"]
                        and interval_1["start"] <= interval_2["end"]
                        and interval_1["end"] >= interval_2["start"]
                ):
                    # Calculate the intersection of intervals
                    start = max(interval_1["start"], interval_2["start"])
                    end = min(interval_1["end"], interval_2["end"])
                    temp_shared_intervals.append({"start": start,
                                                  "end": end,
                                                  "haplotype": interval_1["haplotype"],
                                                  "chromosome": interval_1["chromosome"],
                                                  "certainty_level": 1})
        # Update shared_intervals with the current shared intervals
        shared_intervals = temp_shared_intervals
    return shared_intervals


def non_shared_intervals(interval_lists):
    """
    This function will add intervals with opposite haplotypes, adding a field
    called certainty_level to be -1
    """
    result = []
    interval_dict = defaultdict(list)
    for child_intervals in interval_lists:
        for interval in child_intervals:
            interval_dict[interval["chromosome"]].append(interval)
    for chromosome, intervals in interval_dict.items():
        for i in range(len(intervals) - 1):
            interval1 = intervals[i]
            interval2 = intervals[i + 1]
            start_position = max(interval1["start"], interval2["start"])
            end_position = min(interval1["end"], interval2["end"])
            if start_position <= end_position and interval1["haplotype"] != interval2["haplotype"]:
                new_interval = {
                    "start": start_position,
                    "end": end_position,
                    "chromosome": chromosome,
                    "certainty_level": -1
                }
                result.append(new_interval)
    return result


def plot_interval(interval_list, plot_title, save_dir):
    """
    This function plots intervals as straight lines using Matplotlib.
    """
    fig, ax = plt.subplots()
    for i, interval in enumerate(interval_list):
        start_position = interval["start"]
        end_position = interval["end"]
        haplotype = interval["haplotype"]
        # Plot lines for each interval
        ax.plot([start_position, end_position], [haplotype, haplotype], label=f'Interval {i + 1}')
    ax.set(xlabel='Chromosome Position', ylabel='Haplotype', title=plot_title)
    # Save the plot
    save_path = os.path.join(save_dir, f'{plot_title.replace(" ", "_")}_plot.png')
    # Create the directory if it doesn't exist
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(save_path)
