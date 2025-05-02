import pandas as pd

filenames = ["bacteria_250205_075um_03", "contour_250205_075um_03"]

for filename in filenames:
    # Replace with your actual file paths
    file1 = f"{filename}_G.csv"
    file2 = f"{filename}_R.csv"

    # Read the two CSV files
    df1 = pd.read_csv(file1) #G
    df2 = pd.read_csv(file2) #R

    if "contour" in filename:
        df1 = df1[["NAME",	"NAME.id", "COORD",	"COORD.x", "COORD.y", "POSITION"]]
        df2 = df2[["NAME",	"NAME.id", "COORD",	"COORD.x", "COORD.y", "POSITION"]]
    else:
        df1 = df1[["NAME",	"NAME.id", "LOCATION.x", "LOCATION.y", "POSITION",	"SHAPE.area", "SHAPE.length", "SHAPE.orientation",	"SHAPE.width"]]
        df2 = df2[["NAME",	"NAME.id", "LOCATION.x", "LOCATION.y", "POSITION",	"SHAPE.area", "SHAPE.length", "SHAPE.orientation",	"SHAPE.width"]]


    df1["channel"] = "G"
    df1 = df1.rename(columns={'SHAPE.width.mean': 'SHAPE.width'})
    #print(df1.columns.tolist())
    df2["channel"] = "R"
    #print(df2.columns.tolist())

    # Stack them one below the other
    combined_df = pd.concat([df1, df2], axis=0, ignore_index=True)

    # Optionally save the result
    combined_df.to_csv(f"{filename}_combined.csv", index=False)

print("Files successfully combined and saved to 'combined.csv'")
