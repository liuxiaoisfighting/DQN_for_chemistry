import os
import time

# python CalculateCost v1.0
# by wangjiaze on 2020-Jan-13




def CalculateCost( coordinates,spaceGroup='IMMA', cellParameter=[9.4578, 5.1817,  8.8018, 90, 90, 90]  , mindistance=0, maxdistance=3.7, atomArguments=[4, 0, 0, 25],
                  bondArguments=[3.1, 0, 0, 10], d1_3Arguments=[5, 0, 0, 5]):
    # spaceGroup : string (eg: 'P1', 'IBCA')
    # cellParameter : list (eg: [10, 10, 10, 90, 90, 90]
    # coordinates : Two-dimensional list(eg:[[0.5,0,0], [0,0.5,0],[0,0,0.5],[0.5,0.5,0],[0.5,0,0.5]])

    import CoordinateTransformation as CoorTrans

    def DealCoordinates(coordinates):
        def DealCoordinate(coordinate):
            def Return01(x):
                while x > 1:
                    x = x - 1
                while x < 0:
                    x = x + 1
                return x

            x = coordinate[0]
            y = coordinate[1]
            z = coordinate[2]

            x = Return01(x)
            y = Return01(y)
            z = Return01(z)
            return [x, y, z]

        for i in range(len(coordinates)):
            coordinates[i] = DealCoordinate(coordinates[i])
        return coordinates

    def CalculateLimit27CellCartesian(coordinate, cellParameter, limitDistance):
        x = coordinate[0]
        y = coordinate[1]
        z = coordinate[2]

        # print(limitDistance)
        limit27CellCoordinates = []
        limit27CellCoordinates.append([x, y, z])
        if x <= limitDistance:
            limit27CellCoordinates.append([x+1, y, z])
            if y <= limitDistance:
                limit27CellCoordinates.append([x+1, y+1, z])
                if z <= limitDistance:
                    limit27CellCoordinates.append([x+1, y+1, z+1])
                elif z >= 1-limitDistance:
                    limit27CellCoordinates.append([x+1, y+1, z-1])
            elif y >= 1-limitDistance:
                limit27CellCoordinates.append([x+1, y-1, z])
                if z <= limitDistance:
                    limit27CellCoordinates.append([x+1, y-1, z+1])
                elif z >= 1-limitDistance:
                    limit27CellCoordinates.append([x+1, y-1, z-1])
            if z <= limitDistance:
                limit27CellCoordinates.append([x+1, y, z+1])
            elif z >= 1-limitDistance:
                limit27CellCoordinates.append([x+1, y, z-1])
        elif x >= 1-limitDistance:
            limit27CellCoordinates.append([x-1, y, z])
            if y <= limitDistance:
                limit27CellCoordinates.append([x-1, y+1, z])
                if z <= limitDistance:
                    limit27CellCoordinates.append([x-1, y+1, z+1])
                elif z >= 1-limitDistance:
                    limit27CellCoordinates.append([x-1, y+1, z-1])
            elif y >= 1-limitDistance:
                limit27CellCoordinates.append([x-1, y-1, z])
                if z <= limitDistance:
                    limit27CellCoordinates.append([x-1, y-1, z+1])
                elif z >= 1-limitDistance:
                    limit27CellCoordinates.append([x-1, y-1, z-1])
            if z <= limitDistance:
                limit27CellCoordinates.append([x-1, y, z+1])
            elif z >= 1-limitDistance:
                limit27CellCoordinates.append([x-1, y, z-1])
                
        if y <= limitDistance:
            limit27CellCoordinates.append([x, y+1, z])
            if z <= limitDistance:
                limit27CellCoordinates.append([x, y+1, z+1])
            elif z >= 1-limitDistance:
                limit27CellCoordinates.append([x, y+1, z-1])
        elif y >= 1-limitDistance:
            limit27CellCoordinates.append([x, y-1, z])
            if z <= limitDistance:
                limit27CellCoordinates.append([x, y-1, z+1])
            elif z >= 1-limitDistance:
                limit27CellCoordinates.append([x, y-1, z-1])

        if z <= limitDistance:
            limit27CellCoordinates.append([x, y, z+1])
        elif z >= 1-limitDistance:
            limit27CellCoordinates.append([x, y, z-1])
        #print(limit27CellCoordinates)
        return limit27CellCoordinates

    def CalculateBonds(order_allCoordinate, allLimit27Coordinate_order, cellParameter, mindistance2, maxdistance2):
        cartesian_connectCartesians = {}
        order_cartesian = {}
        allLimit27Cartesian_order = {}

        for order, coordinate in order_allCoordinate.items():
            order_cartesian[order] = CoorTrans.Fractional2Cartesian(coordinate, cellParameter)
        for coordinate, order in allLimit27Coordinate_order.items():
            allLimit27Cartesian_order[tuple(CoorTrans.Fractional2Cartesian(coordinate, cellParameter))] = order
        #print(allLimit27Cartesian_order)
        #print(order_cartesian)
        for order1, cartesian1 in order_cartesian.items():
            cartesian_connectCartesians[tuple(cartesian1)] = []
            for cartesian2, order2 in allLimit27Cartesian_order.items():
                if cartesian1 != tuple(cartesian2):
                    distance2 = CoorTrans.CalculateDistance2(cartesian1, cartesian2)
                    if distance2 > 0 and distance2 >= mindistance2 and distance2 <= maxdistance2:
                        cartesian_connectCartesians[tuple(cartesian1)].append(cartesian2)

        return cartesian_connectCartesians

    def Cost(cartesian_connectNumber, cartesian_connectBonds, cartesian_connectBonds13, cartesian_multiplicity, atomArguments, bondArguments, d1_3Arguments):
        def CalculateAtomPart(cartesian_connectNumber, atomArguments, cartesian_multiplicity):
            atomPartCost = 0
            a = atomArguments[0]
            b = atomArguments[1]
            c = atomArguments[2]
            weight = atomArguments[3]
            sumMultiplicities = 0
            for cartesian, connectNumber in cartesian_connectNumber.items():
                #print(cartesian)
                #print(connectNumber)
                multiplicity = cartesian_multiplicity[cartesian]
                atomPartCost = atomPartCost + weight*(connectNumber-a)**2 * multiplicity
                sumMultiplicities += multiplicity

            if len(cartesian_connectNumber) == 0:
                atomPartCost = 0
            else:
                atomPartCost = atomPartCost/sumMultiplicities
            return atomPartCost

        def CalculateBondPart(cartesian_connectBonds, bondArguments, cartesian_multiplicity):
            bondAllCost = 0
            a = bondArguments[0]
            b = bondArguments[1]
            c = bondArguments[2]
            weight = bondArguments[3]
            bondNumber = 0
            sumMultiplicities = 0
            atomNumber = 0
            for cartesian, connectBonds in cartesian_connectBonds.items():
                atomCost = 0
                multiplicity = cartesian_multiplicity[cartesian]
                atomNumber += multiplicity
                for connectBond in connectBonds:
                    bondCost = (connectBond-a)**2
                    bondNumber += 1
                    sumMultiplicities += multiplicity
                    atomCost += bondCost

                if len(connectBonds) == 0:
                    atomCost = 0
                else:
                    atomCost = atomCost/len(connectBonds)
                bondAllCost += atomCost*multiplicity

            if bondNumber == 0:
                bondPartCost = 0
            else:
                bondPartCost = weight*bondAllCost/atomNumber
            return bondPartCost

        def Calculated1_3Part(cartesian_connectBonds13, d1_3Arguments, cartesian_multiplicity):
            d1_3AllCost = 0
            a = d1_3Arguments[0]
            b = d1_3Arguments[1]
            c = d1_3Arguments[2]
            weight = d1_3Arguments[3]
            bondNumber = 0
            sumMultiplicities = 0
            atomNumber = 0
            for cartesian, connectBonds13 in cartesian_connectBonds13.items():
                atomCost = 0
                multiplicity = cartesian_multiplicity[cartesian]
                atomNumber += multiplicity
                for connectBond13 in connectBonds13:
                    d1_3Cost = (connectBond13-a)**2
                    bondNumber += 1
                    sumMultiplicities += multiplicity
                    atomCost += d1_3Cost
                if len(connectBonds13) == 0:
                    atomCost =0
                else:
                    atomCost = atomCost / len(connectBonds13)
                d1_3AllCost += atomCost * multiplicity

            if bondNumber == 0:
                d1_3PartCost = 0
            else:
                d1_3PartCost = weight * d1_3AllCost / atomNumber
            return d1_3PartCost
        part_cost=[]
        atomPart = CalculateAtomPart(cartesian_connectNumber, atomArguments, cartesian_multiplicity)
        bondPart = CalculateBondPart(cartesian_connectBonds, bondArguments, cartesian_multiplicity)
        d1_3Part = Calculated1_3Part(cartesian_connectBonds13, d1_3Arguments, cartesian_multiplicity)
        part_cost.append(atomPart)
        part_cost.append(bondPart)
        part_cost.append(d1_3Part)
        # print('atomPart: %f' %(atomPart))
        # print('bondPart: %f' %(bondPart))
        # print('d1_3Part: %f' %(d1_3Part))
        cost = round(atomPart + bondPart + d1_3Part, 6)
        part_cost.append(cost)
        return part_cost

    #print(coordinates)
    mindistance2 = mindistance*mindistance
    maxdistance2 = maxdistance*maxdistance

    coordinates = DealCoordinates(coordinates)


    cartesian_multiplicity = {}
    matrices = CoorTrans.GetGeneralMatrix(spaceGroup)
    #print(matrices)
    for coordinate in coordinates:
        cartesian_multiplicity[tuple(CoorTrans.Fractional2Cartesian(coordinate, cellParameter))] = len(CoorTrans.SymmetricOperation(coordinate, matrices))
    allCoordinates = CoorTrans.UniqueAtom2AllAtom(coordinates, spaceGroup)
    #print(allCoordinates)
    order_allCoordinate = {}
    allCoordinate_order = {}
    order_coordinate = {}

    allLimit27Coordinate_order = {}
    i = 1
    for coordinate in coordinates:
        order_coordinate[i] = coordinate
        i += 1
    i = 1
    for coordinate in allCoordinates:
        allCoordinate_order[tuple(coordinate)] = i
        order_allCoordinate[i] = coordinate
        limit27CellCoordinates =  CalculateLimit27CellCartesian(coordinate, cellParameter, maxdistance/ min(cellParameter[0], cellParameter[1], cellParameter[2]))
        for limit27CellCoordinate in limit27CellCoordinates:
            allLimit27Coordinate_order[tuple(limit27CellCoordinate)] = i
        i += 1
    #print(allLimit27Coordinate_order )
    cartesian_connectCartesians = CalculateBonds(order_coordinate, allLimit27Coordinate_order, cellParameter, mindistance2, maxdistance2)
    #print(cartesian_connectCartesians )
    cartesian_connectNumber = {}
    cartesian_connectBonds = {}
    cartesian_connectBonds13 = {}

    for cartesian, connectCartesians in cartesian_connectCartesians.items():
        cartesian_connectNumber[cartesian] = len(connectCartesians)
        connectBonds = []
        connectBonds13 = []
        for connectCartesian in connectCartesians:
            connectBonds.append(CoorTrans.CalculateDistance(cartesian, connectCartesian))
        cartesian_connectBonds[cartesian] = connectBonds

        for i in range(len(connectCartesians)):
            for j in range(i+1, len(connectCartesians)):
                connectBonds13.append(CoorTrans.CalculateDistance(connectCartesians[i], connectCartesians[j]))
        cartesian_connectBonds13[cartesian] = connectBonds13

    cost = Cost(cartesian_connectNumber, cartesian_connectBonds, cartesian_connectBonds13, cartesian_multiplicity, atomArguments, bondArguments, d1_3Arguments)
    print(cost[3])
    return cost


if __name__ == '__main__':
    # print()
    # print('################################')
    # print('              CalculateCost v1.0')
    # print('                     2020-Jan-13')
    # print('                    by wangjiaze')
    # print('################################')
    # print()

    pass
   # CalculateCost('R-3M', [10,10,10,90,90,90], [[0.7779502723802081,0.288100716049866,0.7813547336294311], [0.9662539506822094,0.9965234515119386,0.23910643632747064],[0.4451812617337756,0.1302219280116259,0.23447408750505638]])
