import markdown
from Bio import SeqIO
import regex, re
from Bio.Seq import Seq
from collections import defaultdict
import numpy as np
import pandas as pd
from io import BytesIO
import base64
import tabulate


class SangerBaseCall:
    def __init__(self, sanger_file):
        self.sanger_data = SeqIO.read(sanger_file, "abi")
        self.sanger_seq = str(self.sanger_data.seq).upper()

        self.labels = {}
        signal_length = self.sanger_data.annotations["abif_raw"]["PLOC1"][-1]
        base_locations = self.sanger_data.annotations["abif_raw"]["PLOC1"]
        self.annotation = {}
        self.location_base_annotation = {}

        # 列出每个碱基位置上分别是什么碱基和他对应的测序峰值位置。
        for i in range(len(self.sanger_seq)):
            self.annotation[i] = [self.sanger_seq[i], base_locations[i]]
            self.location_base_annotation[base_locations[i]] = self.sanger_seq[i]

        # 列出每个信号列中的信号值
        channels = {"G": "DATA9", "A": "DATA10", "T": "DATA11", "C": "DATA12"}  # GATC
        self.trace = {}
        for c in channels:
            self.trace[c] = self.sanger_data.annotations["abif_raw"][channels[c]]

        # 根据信号峰值填入碱基位置
        for i in range(signal_length):
            if i in base_locations:
                self.labels[i] = self.location_base_annotation[i]
            else:
                self.labels[i] = ""

        # 整体的噪音
        self.noise = {"G": [], "A": [], "T": [], "C": []}

    def locTartet(self, target_seq, max_mismatch=2):
        # Find the match site and location
        match_seqs = []
        for i in range(int(max_mismatch)):
            target_is_reversed = False
            target_seq = Seq(target_seq.upper())
            match_seqs = regex.findall("(" + str(target_seq) + ")" + '{e<=' + str(i) + '}', self.sanger_seq)

            if match_seqs == []:
                # target_seq = target_seq.reverse_complement()
                match_seqs = regex.findall("(" + str(target_seq.reverse_complement()) + ")" + '{e<=' + str(i) + '}',
                                           self.sanger_seq)
                target_is_reversed = True
            else:
                break

            if match_seqs == []:
                continue
            else:
                break

        if match_seqs == []:
            return

        match_seq = match_seqs[0]
        span = re.search(str(match_seq), self.sanger_seq)
        if span:
            span = span.span()
            if target_is_reversed:
                location_start = span[1]
                location_end = span[0]
            else:
                location_start = span[0]
                location_end = span[1]
            # print(i)
            return (str(match_seq), location_start, location_end)

        else:
            return

    def plotTarget(self, target_seq, plot_window=10):
        from matplotlib import pyplot as plt
        plt.cla()
        plt.close("all")
        # plt.figure(dpi=150)
        # plt.figure(figsize=(20, 4))
        plt.rcParams['figure.figsize'] = (15, 4)
        colors = ("black", "green", "red", "blue")  # GATC
        target_location = self.locTartet(target_seq)
        print(target_location)

        if target_location[1] < target_location[2]:
            target_reverse = False
            base_start = target_location[1] - plot_window
            base_end = target_location[2] + plot_window
        else:
            target_reverse = True
            base_start = target_location[2] - plot_window
            base_end = target_location[1] + plot_window

        if base_start * base_end <= 0:
            print("Plot window size out of range")

        plot_start = self.annotation[base_start][1]
        plot_end = self.annotation[base_end][1]

        # GATC trace
        # self.trace["G"]
        trace_G = self.trace["G"][plot_start:plot_end]
        trace_A = self.trace["A"][plot_start:plot_end]
        trace_T = self.trace["T"][plot_start:plot_end]
        trace_C = self.trace["C"][plot_start:plot_end]

        labels = list(self.labels.values())[plot_start: plot_end]

        fig, ax = plt.subplots()

        for base_i in range(len(labels)):  # 在图上标注碱基
            base = labels[base_i]
            if base != "":
                plt.text(base_i - 1, -100, base)

        plt.tick_params(axis="x",
                        length=0,  # tick长度
                        width=0,  # tick的宽度
                        labelsize=0,  # 刻度值大小
                        bottom=False

                        )
        # plt.xticks(ticks=range(plot_end - plot_start), labels=labels, label="sequence")
        plt.xlabel("Sequence")
        plt.ylabel('Sanger Signal Intensity')
        G = ax.plot(trace_G, color=colors[0], label='G')
        A = ax.plot(trace_A, color=colors[1], label='A')
        T = ax.plot(trace_T, color=colors[2], label='T')
        C = ax.plot(trace_C, color=colors[3], label='C')

        sg_plot_start = self.annotation[target_location[1]][1]
        sg_plot_end = self.annotation[target_location[2]][1]

        if sg_plot_start < sg_plot_end:
            sg_plot_end = self.annotation[target_location[2] - 1][1]
            SG = ax.plot([sg_plot_start - plot_start, sg_plot_end - plot_start], [-150, -150], color="brown",
                         label='sgRNA', marker=">")



        else:
            sg_plot_start = self.annotation[target_location[1]][1]
            SG = ax.plot([sg_plot_end - plot_start, sg_plot_start - plot_start - 9], [-150, -150], color="purple",
                         label='sgRNA', marker="<")

        plt.legend()
        return plt

    def getTargetData(self, target_seq, annalyse_window=10):
        """获取靶位点上各个碱基的信号占比，以及靶位点外的噪音平均值"""
        # plt.rcParams['figure.figsize'] = (15, 4)
        # colors = ("black", "green", "red", "blue")  # GATC
        target_location = self.locTartet(target_seq)
        print(target_location)

        if target_location[1] < target_location[2]:
            target_reverse = False
            base_start = target_location[1] - annalyse_window
            base_end = target_location[2] + annalyse_window
            base_5p_loc = target_location[1]
            print(base_start, base_end)
        else:
            target_reverse = True
            base_start = target_location[2] - annalyse_window
            base_end = target_location[1] + annalyse_window
            base_5p_loc = target_location[1]
            print(base_start, base_end)

        if base_start * base_end <= 0:
            print("Plot window size out of range")

        plot_start = self.annotation[base_start][1]
        plot_end = self.annotation[base_end][1]

        # 噪声统计区间
        noise_start_1 = 100
        noise_end_1 = base_start - 20
        noise_start_2 = base_end + 20

        try:
            noise_end_2 = 700
        except Exception as e:
            print("sanger file is too short!")
            noise_end_2 = base_end + 120

        self.noise = {"G": [], "A": [], "T": [], "C": []}  # 用来记录噪音，后续计算平均值和方差
        print(noise_start_1, noise_end_1)
        for i in range(noise_start_1, noise_end_1):
            if i in self.location_base_annotation:
                peak = i
                peak_base = self.location_base_annotation[peak]
                peak_base_num = list(self.location_base_annotation.keys()).index(peak)
                peak_left_border = peak - 3
                peak_right_border = peak + 3

                signal = {}
                signal["G"] = self.trace["G"][peak_left_border:peak_right_border]
                signal["A"] = self.trace["A"][peak_left_border:peak_right_border]
                signal["T"] = self.trace["T"][peak_left_border:peak_right_border]
                signal["C"] = self.trace["C"][peak_left_border:peak_right_border]

                for base in ["G", "A", "T", "C"]:
                    # print(base+peak_base)
                    if base != peak_base:
                        # print(signal[base])
                        self.noise[base].append(signal[base])

        for i in range(noise_start_2, noise_end_2):
            if i in self.location_base_annotation:
                peak = i
                peak_base = self.location_base_annotation[peak]
                peak_base_num = list(self.location_base_annotation.keys()).index(peak)
                peak_left_border = peak - 3
                peak_right_border = peak + 3

                signal = {}
                signal["G"] = self.trace["G"][peak_left_border:peak_right_border]
                signal["A"] = self.trace["A"][peak_left_border:peak_right_border]
                signal["T"] = self.trace["T"][peak_left_border:peak_right_border]
                signal["C"] = self.trace["C"][peak_left_border:peak_right_border]

                for base in ["G", "A", "T", "C"]:
                    # print(base+peak_base)
                    if base != peak_base:
                        # print(signal[base])
                        self.noise[base].append(signal[base])

        result = {}
        # fig, ax = plt.subplots()
        sheet = pd.DataFrame(index=["G", "A", "T", "C", "sequence", "guide"])
        for i in range(plot_start, plot_end):
            if i in self.location_base_annotation:
                peak = i
                peak_base = self.location_base_annotation[peak]
                peak_base_num = list(self.location_base_annotation.keys()).index(peak)
                peak_trace = self.trace[peak_base]
                # print(peak_trace)eak_val = peak_trace[peak]

                """
                注释掉这一段。这段代码用于查找整个信号的汉宁窗口，然后再以该窗口计算信号、噪音。
                但是，窗口边界可能与其它碱基的窗口重叠，从而导致噪音误判。

                只采用信号峰值附近的应该就够了。
                peak_left_data = np.array(peak_trace[peak-10:peak][::-1])
                peak_right_data = np.array(peak_trace[peak:peak+10])


                peak_left_distance = argrelextrema(peak_left_data,comparator=np.less_equal, order=3)
                peak_right_distance = argrelextrema(peak_right_data,comparator=np.less_equal, order=3)
                peak_left_width = (peak_left_distance[0][0])
                peak_right_width = (peak_right_distance[0][0])
                peak_left_border = peak - peak_left_width
                peak_right_border = peak + peak_right_width

                """

                peak_left_border = peak - 3
                peak_right_border = peak + 3

                signal_G = self.trace["G"][peak_left_border:peak_right_border]
                signal_A = self.trace["A"][peak_left_border:peak_right_border]
                signal_T = self.trace["T"][peak_left_border:peak_right_border]
                signal_C = self.trace["C"][peak_left_border:peak_right_border]

                signal_intensity = {}
                signal_intensity["G"] = sum(signal_G)
                signal_intensity["A"] = sum(signal_A)
                signal_intensity["T"] = sum(signal_T)
                signal_intensity["C"] = sum(signal_C)
                signal_intensity["ALL"] = sum(
                    [signal_intensity["G"], signal_intensity["A"], signal_intensity["T"], signal_intensity["C"]])

                if target_reverse:
                    # print("base_5p_loc",base_5p_loc)
                    peak_base_relative_num = base_5p_loc - peak_base_num

                    if peak_base_relative_num <= 0:
                        peak_base_relative_num = base_5p_loc - peak_base_num - 1
                else:
                    peak_base_relative_num = peak_base_num - base_5p_loc
                    if peak_base_relative_num >= 0:
                        peak_base_relative_num = peak_base_relative_num + 1

                result[peak_base_relative_num] = [peak_base, peak_base_num, signal_intensity]

        # plot_sequence = pd.DataFrame()
        for location in result:
            base = result[location][0]
            sheet.loc["G", location] = (result[location][2]["G"] / result[location][2]["ALL"] * 100)
            sheet.loc["A", location] = (result[location][2]["A"] / result[location][2]["ALL"] * 100)
            sheet.loc["T", location] = (result[location][2]["T"] / result[location][2]["ALL"] * 100)
            sheet.loc["C", location] = (result[location][2]["C"] / result[location][2]["ALL"] * 100)
            sheet.loc["sequence", location] = base

            if (location < len(target_seq)) and (location > 0) and (target_reverse == False):
                sheet.loc["guide", location] = target_seq[location - 1]
            elif (location < len(target_seq)) and (location > 0) and (target_reverse == True):
                sheet.loc["guide", location] = target_seq[location - 1]

            else:
                sheet.loc["guide", location] = ""

        return sheet

    def plotTargetHeatmap(self, data_sheet, base_of_interest="", color_map="YlGn"):
        from matplotlib import pyplot as plt
        plt.cla()
        plt.close("all")
        sequence = data_sheet.loc["sequence"]
        columns = ["G", "A", "T", "C", "Guide"]
        sg_data = []

        # 给sgRNA的位置也赋值，从而带上颜色
        for base in data_sheet.loc["guide"]:
            if base == "":
                sg_data.append(-15)
            else:
                sg_data.append(-50)

        data = np.array([data_sheet.loc["G"],
                         data_sheet.loc["A"],
                         data_sheet.loc["T"],
                         data_sheet.loc["C"],
                         sg_data]).astype("float")

        fig, ax = plt.subplots(figsize=(20, 6))
        ax.set_title(
            "Each base (ACGT) signal intensity at each position along " + "".join(list(data_sheet.loc["guide"])))
        im = ax.imshow(data, cmap=color_map)

        ax.set_xticks(np.arange(len(sequence)), labels=sequence)
        ax.set_yticks(np.arange(len(columns)), labels=columns)

        ax.tick_params(top=True, bottom=True, labeltop=True, labelbottom=True)

        for i in range(len(columns)):
            for j in range(len(sequence)):
                if i < 4:
                    text = ax.text(j, i, round(data[i, j], 1), ha="center", va="center", color="black")
                else:
                    text = ax.text(j, i, list(sequence.index)[j], ha="center", va="center", color="black")
        if base_of_interest == "":
            pass
        elif base_of_interest == "G":
            y = 0
            ax.add_patch(plt.Rectangle((-0.5, y - 0.5), len(sequence), 1, fill=False, edgecolor="red", lw=5))
        elif base_of_interest == "A":
            y = 1
            ax.add_patch(plt.Rectangle((-0.5, y - 0.5), len(sequence), 1, fill=False, edgecolor="red", lw=5))
        elif base_of_interest == "T":
            y = 2
            ax.add_patch(plt.Rectangle((-0.5, y - 0.5), len(sequence), 1, fill=False, edgecolor="red", lw=5))
        elif base_of_interest == "C":
            y = 3
            ax.add_patch(plt.Rectangle((-0.5, y - 0.5), len(sequence), 1, fill=False, edgecolor="red", lw=5))

        return plt

    def plt2md(self, plt):
        buffer = BytesIO()
        plt.savefig(buffer, format='png')

        # 转换base64并以utf8格式输出
        pic_base64 = base64.b64encode(buffer.getvalue()).decode('utf8')
        md = """![](data:image/png;base64,""" + pic_base64 + """)"""
        return md

    def generateReport(self, target_seq, report_name, save_path, plot_window=5, base_of_interest=""):
        sanger_plot = self.plotTarget(target_seq=target_seq, plot_window=plot_window)
        sanger_plot_md = self.plt2md(sanger_plot)

        perc_around_guide_sheet = self.getTargetData(target_seq=target_seq, annalyse_window=plot_window)
        heatmap_plot = self.plotTargetHeatmap(perc_around_guide_sheet, base_of_interest)
        heatmap_plot_md = self.plt2md(heatmap_plot)

        title = "# BaseEdit Report of " + report_name + "\n"
        sep = "\n\n---------"
        description = """## Description 
        \n\nthis file contains the sanger sequencing report around the guide 
        """ + target_seq + "\n"

        sub_title_1 = "\n\n## Sanger sequencing result around " + target_seq + "\n\n"
        Fig1_annotation = "\n\nThis figure shows the Sanger sequencing signal quality."

        sub_title_2 = "\n\n## Base call result in each position around " + target_seq + "\n"
        Fig2_annotation = """\n\nThis figure shows the possibilities of each base at each sequencing position. 
                                Guide position at the botton line shows the relative position to your guide RNA.\n"""
        if base_of_interest != "":
            Fig2_annotation = Fig2_annotation + "\n\nThe red rectangle is the base you want to convert to."
            description = description + "\n\nThis sample is supposed to be converted to " + base_of_interest

        sub_title_3 = "\n\n## Base call table in each position around " + target_seq + "\n"
        sheet = perc_around_guide_sheet.to_markdown()
        sheet_annotation = """\n\n This sheet shows the possibilities of each base at each sequencing position.\n\nIt's the same as the heatmap above.
                                \n\nYou can copy the entire sheet and pasted it to Excel for further analysis.\n\n"""
        content = [title, description, sep, sub_title_1, sanger_plot_md, Fig1_annotation, sep, sub_title_2,
                   heatmap_plot_md, Fig2_annotation, sep, sub_title_3, sheet_annotation, sheet]

        html_content = markdown.markdown("".join(content), extensions=["tables"])
        with open(save_path, "w") as f:
            f.write(html_content)


if __name__ == "__main__":

    a = SangerBaseCall("example.ab1")
    a.generateReport("CACTGGAATGACACACGCCC", "report_name", "test.html",plot_window=5)
