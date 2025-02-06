import os
import re
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main(directory_path):
    # Initialize lists to store data
    first_ligand_scores = []
    lowest_top_5_ligand_scores = []
    lowest_file_rmsd_scores = []
    below_2_5_count_first = 0
    below_2_5_count_lowest_top_5 = 0
    below_2_5_count_lowest_file = 0

    # Iterate through each folder in the directory
    for folder in os.listdir(directory_path):
        folder_path = os.path.join(directory_path, folder)

        # Check if it's a directory
        if os.path.isdir(folder_path):
            # Find the corresponding rmsd_output.txt file
            file_path = os.path.join(folder_path, f"{folder}_rmsd_output.txt")

            # Check if the file exists
            if os.path.exists(file_path):
                # Read the file and extract the first ligand score, lowest among top 5 ligands, and lowest from the file
                with open(file_path, 'r') as file:
                    content = file.read()
                    matches = re.findall(r'RMSD for molecule (\d+): (\d+\.\d+)', content)

                    if matches:
                        # Extract scores for the first ligand, the top 5 ligands, and the lowest from the file
                        first_ligand_score = float(matches[0][1])
                        lowest_top_5_ligand_score = min(float(match[1]) for match in matches[:5])
                        lowest_file_rmsd_score = min(float(match[1]) for match in matches)

                        first_ligand_scores.append(first_ligand_score)
                        lowest_top_5_ligand_scores.append(lowest_top_5_ligand_score)
                        lowest_file_rmsd_scores.append(lowest_file_rmsd_score)

                        # Check if the scores are below 2.0
                        if first_ligand_score < 2.0:
                            below_2_5_count_first += 1

                        if lowest_top_5_ligand_score < 2.0:
                            below_2_5_count_lowest_top_5 += 1

                        if lowest_file_rmsd_score < 2.0:
                            below_2_5_count_lowest_file += 1

    # Calculate the percentages of scores below 2.0
    percentage_below_2_5_first = (below_2_5_count_first / len(first_ligand_scores)) * 100 if first_ligand_scores else 0
    percentage_below_2_5_lowest_top_5 = (below_2_5_count_lowest_top_5 / len(lowest_top_5_ligand_scores)) * 100 if lowest_top_5_ligand_scores else 0
    percentage_below_2_5_lowest_file = (below_2_5_count_lowest_file / len(lowest_file_rmsd_scores)) * 100 if lowest_file_rmsd_scores else 0

    # Create pandas DataFrames
    data_first = pd.DataFrame({'First Ligand Score': first_ligand_scores})
    data_lowest_top_5 = pd.DataFrame({'Lowest Top 5 Ligand Score': lowest_top_5_ligand_scores})
    data_lowest_file = pd.DataFrame({'Lowest File RMSD Score': lowest_file_rmsd_scores})

    # Generate a combined violin plot
    plt.figure(figsize=(12, 6))
    sns.violinplot(data=[data_first['First Ligand Score'], data_lowest_top_5['Lowest Top 5 Ligand Score'],
                         data_lowest_file['Lowest File RMSD Score']],
                   palette=['skyblue', 'orange', 'green'])
    plt.ylabel('RMSD relative to native pose', fontsize=14)
    plt.xticks([0, 1, 2], ['Top ranked pose', 'Lowest RMSD pose of top 5 ', 'Lowest RMSD pose of 100'], fontsize=14)
    plt.ylim(-2, 12)  # Set the limits of the y-axis to go from 0 to 12
    plt.show()

    # Display results for first ligand scores
    if first_ligand_scores:
        min_score_first = min(first_ligand_scores)
        max_score_first = max(first_ligand_scores)
        print(f"Minimum First Ligand Score: {min_score_first:.2f}")
        print(f"Maximum First Ligand Score: {max_score_first:.2f}")
        print(f"Percentage of First Ligand scores below 2.0: {percentage_below_2_5_first:.2f}%")

    # Display results for lowest among top 5 ligands
    if lowest_top_5_ligand_scores:
        min_score_lowest_top_5 = min(lowest_top_5_ligand_scores)
        max_score_lowest_top_5 = max(lowest_top_5_ligand_scores)
        print(f"\nMinimum Lowest Top 5 Ligand Score: {min_score_lowest_top_5:.2f}")
        print(f"Maximum Lowest Top 5 Ligand Score: {max_score_lowest_top_5:.2f}")
        print(f"Percentage of Lowest Top 5 Ligand scores below 2.0: {percentage_below_2_5_lowest_top_5:.2f}%")

    # Display results for lowest file RMSD scores
    if lowest_file_rmsd_scores:
        min_score_lowest_file = min(lowest_file_rmsd_scores)
        max_score_lowest_file = max(lowest_file_rmsd_scores)
        print(f"\nMinimum Lowest File RMSD Score: {min_score_lowest_file:.2f}")
        print(f"Maximum Lowest File RMSD Score: {max_score_lowest_file:.2f}")
        print(f"Percentage of Lowest File RMSD scores below 2.0: {percentage_below_2_5_lowest_file:.2f}%")

if __name__ == "__main__":
    # Command-line argument for directory path
    if len(sys.argv) != 2:
        print("Usage: python rmsd_trends.py <parent_directory>")
        sys.exit(1)

    input_directory = sys.argv[1]

    # Validate the provided directory
    if not os.path.isdir(input_directory):
        print(f"Error: {input_directory} is not a valid directory.")
        sys.exit(1)

    # Run the main function
    main(input_directory)
