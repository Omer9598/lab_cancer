import os

import matplotlib.pyplot as plt


def create_intervals(haplotype_dict: dict, interval_len=1000000):
    """
    This function will create a list for each child in the following format:
    [(interval number is the index) {start position: , end position: ,
    haplotype: (1 or 2)}
    each interval starts with the position of a variant from haplotype 1 or 2,
    and ends when the next variant is from the opposite haplotype, where a new
    interval will start
    """
    positions = sorted(haplotype_dict.keys())
    len_positions = len(positions)
    intervals = []
    index = -1
    while index < len_positions - 1:
        index += 1
        interval_start = positions[index]
        interval_end = None
        cur_position = interval_start
        next_position = interval_start
        chromosome = haplotype_dict[cur_position][-3]
        cur_haplotype = haplotype_dict[cur_position][-2]
        next_haplotype = haplotype_dict[next_position][-2]

        while index < len_positions - 1:
            # Distance between 2 variants is too big
            if next_position - cur_position > interval_len:
                interval_end = cur_position + interval_len
                break
            # Different haplotype
            if cur_haplotype != next_haplotype:
                interval_end = cur_position
                break
            # None of the conditions met, extending the interval
            cur_position = positions[index]
            index += 1
            next_position = positions[index]
            next_haplotype = haplotype_dict[next_position][-2]
        # Finished with the current interval - adding to the list
        if not interval_end:
            interval_end = positions[-1]
        current_interval = {"start": interval_start,
                            "end": interval_end,
                            "haplotype": cur_haplotype,
                            "chromosome": chromosome}
        intervals.append(current_interval)
    return intervals


def shared_interval(interval_lists):
    """
    This function creates a new list containing intervals that are shared
    in all the lists given, according to the haplotype.
    It also adds certainty level of 1 to the shared intervals.
    If haplotypes are different, haplotype is set to 0 and certainty level to -1.
    """
    shared_intervals = interval_lists[0]
    new_shared_intervals = []
    for interval_list in interval_lists[1:]:
        for interval_1 in shared_intervals:
            for interval_2 in interval_list:
                # Check if intervals overlap
                if interval_1["start"] <= interval_2["end"] and interval_1["end"] >= interval_2["start"]:
                    cur_shared_interval = {
                        "start": max(interval_1["start"], interval_2["start"]),
                        "end": min(interval_1["end"], interval_2["end"]),
                        "haplotype": interval_1["haplotype"] if interval_1["haplotype"] == interval_2[
                            "haplotype"] else 0,
                        "chromosome": interval_1["chromosome"],
                        "certainty_level": 1 if interval_1["haplotype"] == interval_2["haplotype"] else -1
                    }
                    new_shared_intervals.append(cur_shared_interval)
        shared_intervals = new_shared_intervals
        new_shared_intervals = []

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
        color = 'red' if certainty_level == -1 else 'green'
        ax.plot([start_position, end_position], [haplotype, haplotype], label=f'Interval {i + 1}', color=color)
    ax.set(xlabel='Chromosome Position', ylabel='Haplotype', title=plot_title)
    save_path = os.path.join(save_dir, f'{plot_title.replace(" ", "_")}_plot.png')
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(save_path)
