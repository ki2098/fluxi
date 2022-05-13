def build_cell(axis):
    f = open(axis+".txt")
    line = f.readline()
    num_grid = int(line)
    grid = [0.0] * num_grid
    
    line = f.readline()
    i = 0
    while line:
        grid[i] = float(line)
        i += 1
        line = f.readline()
    f.close()
    
    num_cell = num_grid - 1
    cell = [0.0] * (num_cell + 4)
    for i in range(2, num_cell + 2):
        cell[i] = 0.5 * (grid[i - 1] + grid[i - 2])
    
    cell[1] = 2 * grid[0] - cell[2]
    cell[0] = 2 * cell[1] - cell[2]
    cell[num_cell + 2] = 2 * grid[num_grid - 1] - cell[num_cell + 1]
    cell[num_cell + 3] = 2 * cell[num_cell + 2] - cell[num_cell + 1]

    return cell, num_cell

if __name__ == "__main__":
    cell, num_cell = build_cell("x")
    print(cell)

    f = open("x.cell", "w")
    for i in range(num_cell + 4):
        f.write("%s\n"%('%.15E'%(cell[i])))
    f.close

    cell, num_cell = build_cell("y")
    print(cell)

    f = open("y.cell", "w")
    for i in range(num_cell + 4):
        f.write("%s\n"%('%.15E'%(cell[i])))
    f.close

    cell, num_cell = build_cell("z")
    print(cell)

    f = open("z.cell", "w")
    for i in range(num_cell + 4):
        f.write("%s\n"%('%.15E'%(cell[i])))
    f.close