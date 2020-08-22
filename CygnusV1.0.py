# Sagittarius clone

# controls:
# the A and D keys are used for movement
# you click mb1 and drag your mouse *away* from where you want to shoot
# An indicator shows where your arrow is going to land
# let go of mb1 to shoot


# distance between two points
def dist(v1, v2):
	return ((v1[0]-v2[0])**2 + (v1[1]-v2[1])**2)**0.5

# scales vectors by a constant
def scale(v, k):
	for i in range(2):
		v[i] *= k/magnitude(v)
	return v

# returns the magnitude of a vector
def magnitude(v):
	return (v[0]**2+v[1]**2)**0.5

# sign function
def sgn(x):
    return x//abs(x)


# finds the inverse of a 2x2 matrix
def inverseMatrix2(rectMatrix):
	basis1 = rectMatrix[0]
	basis2 = rectMatrix[1]

	matrix = [[basis1[0], basis2[0]], [basis1[1], basis2[1]]]

	determinant = matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]
	if determinant == 0:
		return 'invalid'

	inverse = [[basis2[1], -basis2[0]], [-basis1[1], basis1[0]]]
	for i in range(2):
		for j in range(2):
			inverse[i][j] /= determinant

	return inverse


# finds the product of a 2x2 and 2x1 matrix
def matrixProd(A, V):
	prod = []
	prod.append(A[0][0]*V[0]+A[0][1]*V[1])
	prod.append(A[1][0]*V[0]+A[1][1]*V[1])

	return prod


# checks for intersection
def check(rectVectors, pVector):
	rectMatrix = []
	for i,j in zip(rectVectors[0], rectVectors[1]):
		rectMatrix.append([i,j])

	Ainv = inverseMatrix2(rectMatrix)
	if Ainv == 'invalid':
		return False

	prod = matrixProd(Ainv, pVector)

	if (prod[0] > 0 and prod[1] > 0) and (prod[0] < 1 and prod[1] < 1):
		return True
	else:
		return False


# checks for collisions
def collidePlayer(rectVectors, arrowVectors):
	basisRect = []
	basisArrow = []

	for i in [1,3]:
		basisRect.append([rectVectors[i][0] - rectVectors[0][0], rectVectors[i][1] - rectVectors[0][1]])
		basisArrow.append([arrowVectors[i][0] - arrowVectors[0][0], arrowVectors[i][1] - arrowVectors[0][1]])

	offsetRect = rectVectors[0]
	offsetArrow = arrowVectors[0]

	for i in range(4):
		arrowPoint = [arrowVectors[i][0] - offsetRect[0], arrowVectors[i][1] - offsetRect[1]]
		rectPoint = [rectVectors[i][0] - offsetArrow[0], rectVectors[i][1] - offsetArrow[1]]

		collide = check(basisRect, arrowPoint) or check(basisArrow, rectPoint)
		if collide:
			return True

	return False


# initializes the game
def initializeGame():
	# dimensions and colour of the canvas
    global canvasWidth, canvasHeight, canvasColour
    windowInfo = pygame.display.Info()
    canvasWidth = 1600
    canvasHeight = 900
    canvasWidth = windowInfo.current_w
    canvasHeight = windowInfo.current_h
    canvasColour = (0,0,0)

    # max radius of the planet
    global maxRadius
    maxRadius = 120

    # distance planet's surfaces have to be from the edge
    global edgeBuffer
    edgeBufferConst = 100
    edgeBuffer = maxRadius + edgeBufferConst

    # distance planet's surfaces have to be apart
    global planetBuffer
    planetBuffer = 130

    # gravitational constant
    # the force scaling that players shoot with
    # the min/max force that players can shoot a projectile
    global gravConst, forceConst, minForce, maxForce
    gravConst = 2000
    forceConst = 0.5
    minForce = 20
    maxForce = 200

    # the cnavas object itself
    global canvas
    canvas = pygame.display.set_mode((canvasWidth, canvasHeight), flags = pygame.FULLSCREEN)

    # time/frame management
    global clock, timeConstant, fps, dt
    clock = pygame.time.Clock()
    timeConstant = 0.004
    fps = 60

    # the projectile and whether it currently exists
    global proj, projExist
    projExist = False

    # an array of player objects and desired length of the playerArray
    global playerArray, numPlayers
    playerArray = []
    numPlayers = 5

    # iframes for the player shooting the projectile
    # how many seconds/frames the projectile lasts before despawning
    global iframes, framesLifeTime
    iframes = 45
    lifeTime = 8
    framesLifeTime = fps*lifeTime

    # a list of positions to draw where previous arrows landed
    # acts as a marker for previous shots
    global phantoms
    phantoms = []

    # indents the phantom arrows further into the planet
    global indentConstant
    indentConstant = 0

    # indents the phantom arrows further into the planet
    global phantomBuffer
    phantomBuffer = 12

    # indents the player's base into the planet
    global standingBuffer
    standingBuffer = 2

    # colour of the turn-indicator
    global indicatorColour
    indicatorColour = (255,0,0)

    # the trail of the projectile as well as how many points can be in it
    global projTrail, trailLength
    projTrail = Trail((0,0,0))
    trailLength = 25

    # a dict between player index and player colours
    global playerColours
    playerColours = {}
    playerColours[0] = (255,0,0)
    playerColours[1] = (0,255,0)
    playerColours[2] = (0,0,255)
    playerColours[3] = (255,165,0)
    playerColours[4] =(208,0,208)

    # array that contains all star objects
    global starArray
    starArray = []

    # whether the force vector should be drawn
    global drawForceVector
    drawForceVector = False

    pygame.mouse.set_visible(False)


# initializes planetArray
def initializePlanets():
	# an array containing all Planet objects
	global planetArray
	numPlanets = 5
	planetArray = []

	# info on what the dimmest/brightest RGB values are
	dimmestColour = 80
	brightestColour = 250

	# seed(34)
	# seed(25)
	# seed(20)

	# position needs to be intiailized
	position = []

	for i in range(numPlanets):
		planetSensible = False

		# this loop runs until a planet that meets all the criteria are found
		while not planetSensible:

			# selects the position and size of the planet
			radius = randrange(maxRadius//5, maxRadius)
			xpos = randrange(radius + edgeBuffer, canvasWidth - radius - edgeBuffer)
			ypos = randrange(radius + edgeBuffer, canvasHeight - radius - edgeBuffer)
			position = [xpos, ypos]

			# checks whether all constraints are met
			planetSensible = True

			for planet in planetArray:
				if not planet.farEnough(position, radius):
					planetSensible = False
					break

		# selects a random colour for the planet
		planetColour = [0, 0, 0]
		for i in range(3):
			planetColour[i] = randrange(dimmestColour, brightestColour)
		planetColour = tuple(planetColour)

		planetArray.append(Planet(radius, position, planetColour))


# initializes the starArray
def initializeStars():
	# min distance from star to a star or planet
	minDistFromBody = 15
	numStars = 50

	while len(starArray) < numStars:
		xpos = randrange(0, canvasWidth)
		ypos = randrange(0, canvasHeight)
		position = [xpos, ypos]

		# track if the min distance constraints are met
		# if they are: create the star
		validStar = True
		for planet in planetArray:
			distance = dist(planet.pos, position)
			if distance - planet.r < minDistFromBody:
				validStar = False
				break

		if validStar:
			for star in starArray:
				distance = dist(star.pos, position)
				if distance - star.r < minDistFromBody:
					validStar = False
					break

			if validStar:
				starArray.append(Star(position))


# draws the force vector when a player is shooting
def drawForce(turn, mouseClickPos):
	# scales the vector so it looks visually decent
	scaling = 0.8
	mouseCurrentPos = pygame.mouse.get_pos()
	
	# normalizes direction vector of force
	# scales it accordingly
	forceDir = [mouseCurrentPos[0] - mouseClickPos[0], mouseCurrentPos[1] - mouseClickPos[1]]

	mag = magnitude(forceDir)
	force = forceConst * mag
	if force < minForce:
		force = minForce

	elif force > maxForce:
		force = maxForce

	if mag:
		forceDir = [forceDir[0]/mag, forceDir[1]/mag]
		unitForce = forceDir.copy()
		forceDir = [scaling*force*forceDir[0], scaling*force*forceDir[1]]
	else:
		return

	# initializes the player object that will be used for colour and position calculations
	player = playerArray[turn]

	# finds where the force vector will end
	forceVector = [player.pos[0] - forceDir[0], player.pos[1] - forceDir[1]]

	# gap is the distance between player.pos and the force vector
	# minDistToDraw is the minimum distance to draw the vector as well as the arrow
	gap = 60
	minDistToDraw = 120
	if mag > minDistToDraw:
		forceVector = [int(forceVector[0]), int(forceVector[1])]
		drawForce = [int(forceVector[0]), int(forceVector[1])]

		drawPos = player.pos.copy()
		
		drawPos = [int(drawPos[0] - gap*unitForce[0]), int(drawPos[1] - gap*unitForce[1])]

		# width of the drawn vector
		lineWidth = 5
		pygame.draw.line(canvas, player.colour, forceVector, drawPos, lineWidth)

	drawVectorHead(player, forceVector, unitForce)


def drawVectorHead(player, forceVector, unitForce):
	# unit normal needed for later calculations
	unitNormal = [unitForce[1], -unitForce[0]]

	# draws an arrow indicating the direction of the force of the shot
	# dimensions of the arrow as well as how far it extrudes from the point
	arrowHeight = 1
	arrowWidth = 10
	arrowExtrude = 10

	arrowVertices = [[0,0],[0,0],[0,0]]
	arrowVertices[0] = forceVector.copy()
	arrowVertices[0] = [arrowVertices[0][0] - arrowExtrude*unitForce[0], arrowVertices[0][1] - arrowExtrude*unitForce[1]]

	arrowBase = [forceVector[0] + arrowHeight*unitForce[0], forceVector[1] + arrowHeight*unitForce[1]]
	arrowVertices[1] = [arrowBase[0] - arrowWidth*unitNormal[0]/2, arrowBase[1] - arrowWidth*unitNormal[1]]
	arrowVertices[2] = [arrowBase[0] + arrowWidth*unitNormal[0]/2, arrowBase[1] + arrowWidth*unitNormal[1]]

	# integers work better with pygame.draw
	for i in range(3):
		for j in range(2):
			arrowVertices[i][j] = int(arrowVertices[i][j])

	pygame.draw.polygon(canvas, player.colour, arrowVertices)


# a loop that draws every frame
def drawLoop(turn, mouseClickPos):
	canvas.fill(canvasColour)

	for planet in planetArray:
		planet.draw()

	for star in starArray:
		star.draw()

	if projExist:
		projTrail.draw()
		proj.draw()

	# only draw the turn indicator over the player whose turn it is
	playerArray[turn].drawIndicator()

	for player in playerArray:
		if player.playerExist:
			player.draw()

	for phantom in phantoms:
		phantom.draw()

	if drawForceVector:
		drawForce(turn, mouseClickPos)

	pygame.display.update()


# finds the net graviational force on the projectile
def findPull(proj):

	# reinitializes the acceleration vector
	proj.accel = [0,0]
	rx0 = proj.pos[0]
	ry0 = proj.pos[1]

	# calculates and accumulates it for every planet
	# calculations are based on the graviational force formula
	for planet in planetArray:
		rxp = planet.pos[0]
		ryp = planet.pos[1]

		rx = rxp - rx0
		ry = ryp - ry0

		d = dist(planet.pos, proj.pos)

		# d is cubed in this expression as rx and ry need to be normalize
		# d**1 to normalize, and d**2 due to the inverse square law
		# this gives us d**3 overall
		proj.accel[0] += rx * planet.m / d**3
		proj.accel[1] += ry * planet.m / d**3

	for i in range(2):
		proj.accel[i] *= dt*gravConst


# checks for planet collision
def collidePlanet(proj):
	for i in range(len(planetArray)):

		# simple point-circle intersection algorithm
		if abs(dist(planetArray[i].pos, proj.pos)) < abs(planetArray[i].r):

			# returns whether it has intersected, and with which planet index in planetArray
			return [True, i]
	return [False, False]


class Planet:
	def __init__(self, radius, position, colour):
		self.pos = position
		self.r = radius
		self.m = radius**2

		self.colour = colour

	def draw(self):
		pygame.draw.circle(canvas, self.colour, self.pos, self.r)

	def farEnough(self, position, radius):
		return dist(self.pos, position) > self.r + radius + planetBuffer


class Projectile:
	def __init__(self, position, velocity, colour):
		self.pos = position
		self.v = velocity
		self.accel = [0,0]

		# load the sprite then convert it to a pixel array to change its colour to the player's colour
		# the proj being the player's colour helps them keep track of their shots
		spriteFile = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Sprites\\Arrow20.png")
		self.sprite = pygame.image.load(spriteFile).convert_alpha()
		
		colouredSprite = pygame.PixelArray(self.sprite)
		self.colourDist = 0.52
		self.baseColour = (254,0,0)
		colouredSprite.replace(self.baseColour, colour, distance = self.colourDist)

		self.spriteDim = self.sprite.get_rect().size

		angle = atan2(self.v[1], self.v[0])*180/pi

		# updates all vertice's positions
		angle = atan2(self.v[1], -self.v[0])

		hitScale = 3/4
		x1 = hitScale*self.spriteDim[0]*sin(angle)
		y1 = hitScale*self.spriteDim[0]*cos(angle)

		x2 = y1
		y2 = -x1

		v1 = [self.pos[0] - 0.5*(x1+x2), self.pos[1] - 0.5*(y1+y2)]
		v2 = [v1[0] + x1, v1[1] + y1]
		v3 = [v2[0] + x2, v2[1] + y2]
		v4 = [v1[0] + x2, v1[1] + y2]

		self.vertices = [v1, v2, v3, v4]

		self.colour = colour

	# moves the projectile in accordance to the planets' gravities
	def path(self):
		findPull(self)

		for i in range(2):
			self.v[i] += self.accel[i]*dt
			self.pos[i] += self.v[i]*dt

		# updates all vertice's positions
		angle = atan2(self.v[1], -self.v[0])

		hitScale = 3/4
		x1 = hitScale*self.spriteDim[0]*sin(angle)
		y1 = hitScale*self.spriteDim[0]*cos(angle)

		x2 = y1
		y2 = -x1

		v1 = [self.pos[0] - 0.5*(x1+x2), self.pos[1] - 0.5*(y1+y2)]
		v2 = [v1[0] + x1, v1[1] + y1]
		v3 = [v2[0] + x2, v2[1] + y2]
		v4 = [v1[0] + x2, v1[1] + y2]

		self.vertices = [v1, v2, v3, v4]

		projTrail.updateTrail()

	def draw(self):
		drawPos = []
		for i in range(2):
			drawPos.append(int(self.pos[i] - self.spriteDim[i]/2))

		rotationAngle = atan2(self.v[1], -self.v[0])*180/pi
		newSprite = self.sprite.copy()
		newSprite = pygame.transform.rotate(newSprite, rotationAngle)

		canvas.blit(newSprite, drawPos)

		# debug by drawing projectile vertices
		# pygame.draw.polygon(canvas, (255,255,255), self.vertices)
		# pygame.draw.circle(canvas, [255,255,255], (int(self.vertices[0][0]), int(self.vertices[0][1])), 5)


# the trails that follows projectiles
class Trail:
	def __init__(self, colour):
		# an array of points on the trail
		self.points = []
		self.colour = colour
		self.r = 2

	# adds the projectile's most recent point to the trail then limits its size
	def updateTrail(self):
		newPoint = [int(proj.pos[0]), int(proj.pos[1])]
		self.points.append(newPoint)
		if len(self.points) > trailLength:
			del self.points[0]

	def draw(self):
		# don't draw the most recent 2-3 trails as sometimes they look offset
		for point in self.points[:-2]:
			pygame.draw.circle(canvas, self.colour, point, self.r)


class Player:
	def __init__(self, planetIndex, colour):
		self.angle = pi*randrange(360)/180
		self.planet = planetArray[planetIndex]

		self.walkDist = pi/(self.planet.r*2)
		self.width = 20
		self.height = 60
		self.colour = colour

		self.vertices = [[0,0],[0,0],[0,0],[0,0]]
		self.pos = [0,0]

		self.playerExist = True

		self.indicatorHeight = 10
		self.indicatorWidth = 2

		self.indicatorVerts = [[0,0],[0,0],[0,0]]

		# vertical distance between the player's base and the start of the indicator
		self.indicatorGap = 2*self.height/3

		self.updatePos()

	# updates the position of the player on the planet
	def updatePos(self):
		surfaceVector = [(self.planet.r - standingBuffer) * cos(self.angle), (self.planet.r - standingBuffer) * sin(self.angle)]
		normalVector = [-surfaceVector[1], surfaceVector[0]]

		for i in range(2):
			normalVector[i] *= self.width/(2*self.planet.r)

		for i in range(2):
			self.vertices[0][i] = self.planet.pos[i] + surfaceVector[i] + normalVector[i]
			self.vertices[1][i] = self.planet.pos[i] + surfaceVector[i] - normalVector[i]
			self.vertices[2][i] = self.vertices[1][i] + self.height * surfaceVector[i]/magnitude(surfaceVector) 
			self.vertices[3][i] = self.vertices[0][i] + self.height * surfaceVector[i]/magnitude(surfaceVector)

		for v in range(4):
			for i in range(2):
				self.vertices[v][i] = int(self.vertices[v][i])

		self.pos = [0.5*(self.vertices[0][0]+self.vertices[2][0]),0.5*(self.vertices[0][1]+self.vertices[2][1])]

		self.updateIndicator(surfaceVector, normalVector)

	def updateIndicator(self, vector, unitNormal):
		# normalize the vector
		mag = magnitude(vector)
		for i in range(2):
			vector[i] /= mag

		# vector calculations to find the vertices of an isoceles triangle that is rotated to the player's orientation
		self.indicatorVerts[0][0] = int(vector[0] + (self.indicatorGap) * cos(self.angle) + self.pos[0])
		self.indicatorVerts[0][1] = int(vector[1] + (self.indicatorGap) * sin(self.angle) + self.pos[1])

		self.indicatorVerts[1][0] = self.indicatorVerts[0][0] + vector[0] * self.indicatorHeight + unitNormal[0] * self.indicatorWidth/2
		self.indicatorVerts[1][1] = self.indicatorVerts[0][1] + vector[1] * self.indicatorHeight + unitNormal[1] * self.indicatorWidth/2

		self.indicatorVerts[2][0] = self.indicatorVerts[0][0] + vector[0] * self.indicatorHeight - unitNormal[0] * self.indicatorWidth/2
		self.indicatorVerts[2][1] = self.indicatorVerts[0][1] + vector[1] * self.indicatorHeight - unitNormal[1] * self.indicatorWidth/2
		
	def walk(self, key):
		if key == pygame.K_a:
			self.angle -= self.walkDist
			self.updatePos()
		elif key == pygame.K_d:
			self.angle += self.walkDist
			self.updatePos()

	def shoot(self, origin, pullback):
		global proj, projExist

		netVector = []
		dirVector = []

		for i in range(2):
			netVector.append(origin[i] - pullback[i])
			dirVector.append(netVector[i])

		force = 0 

		# avoids a zero-division error if the player shoots with 0 force
		if magnitude(netVector):
			for i in range(2):
				dirVector[i] /= magnitude(netVector)

			force = magnitude(netVector)

		# constains the force of the players shot to the domain [minForce, maxForce]
		if force < minForce:
			force = minForce
		elif force > maxForce:
			force = maxForce

		resultant = [forceConst*force*dirVector[0], forceConst*force*dirVector[1]]

		# create the projectile object and set its existance to True
		proj = Projectile(self.pos.copy(), resultant, self.colour)

		projExist = True

	def draw(self):
		drawVertices = self.vertices.copy()
		for i in range(4):
			for j in range(2):
				drawVertices[i][j] = int(drawVertices[i][j])

		pygame.draw.polygon(canvas, self.colour, drawVertices)

	def drawIndicator(self):
		drawVertices = self.indicatorVerts.copy()
		for i in range(3):
			for j in range(2):
				drawVertices[i][j] = int(drawVertices[i][j])

		pygame.draw.polygon(canvas, indicatorColour, drawVertices)


# class that keeps track of the arrows previous collisions with planets
# used to draw in arrows when they hit the ground
class PhantomArrow:
	def __init__(self, position, velocity, colour, planetIndex):

		# some vector manipulation that calculates the desirable angle, and position for the phantom arrow relative to the planet
		self.planet = planetArray[planetIndex]
		xcomp = position[0] - self.planet.pos[0]
		ycomp = position[1] - self.planet.pos[1]

		rotationAngle = atan2(velocity[1], - velocity[0])*180/pi
		spriteFile = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Sprites\\Arrow20.png")
		self.sprite = pygame.image.load(spriteFile).convert_alpha()
		
		colouredSprite = pygame.PixelArray(self.sprite)
		colouredSprite.replace(proj.baseColour, colour, distance = proj.colourDist)

		self.sprite = pygame.transform.rotate(self.sprite, rotationAngle)

		normalAngle = atan2(velocity[0], - velocity[1])*180/pi
		positionProjection = [position[0] - self.planet.pos[0], position[1] - self.planet.pos[1]]
		planetAngle = atan2(positionProjection[1], positionProjection[0])

		angleDiff = normalAngle + planetAngle

		hitRadius = self.planet.r - indentConstant*sin(angleDiff)

		self.pos = [int(self.planet.pos[0] + hitRadius*cos(planetAngle) - phantomBuffer*sgn(hitRadius)), int(self.planet.pos[1] + hitRadius*sin(planetAngle)- phantomBuffer*sgn(hitRadius))]

	def draw(self):
		canvas.blit(self.sprite, self.pos)


# visual stars that have no gameplay applications
class Star:
	def __init__(self, position):
		self.pos = position
		smallestR = 2
		largestR = 5
		self.r = randrange(smallestR, largestR)
		self.colour = (250,250,250)

	def draw(self):
		pygame.draw.circle(canvas, self.colour, self.pos, self.r)


# class containing the death particles that are drawn in the wake of a kill
class Particle:
	def __init__(self):
		pass


def mainLoop():
	global dt
	global proj, projExist, projTrail
	global playerColours
	global drawForceVector

	# variables that keep track of whether a player is walking and what direction they are walking in
	walking = False
	walkDir = 0

	# variables that keep track of where the mouse is clicked and released
	mouseClickVector = [0,0]
	mouseReleaseVector = [0,0]

	iframeCount = 0
	framesDuration = 0

	deadPlayers = 0
	playerTurn = randrange(numPlayers)

	mouseClick = False

	for i in range(numPlayers):
		playerArray.append(Player(i, playerColours[i]))

	while True:
		# if everyone except one player is dead, the game is over
		if numPlayers - deadPlayers == 1:
			drawLoop(playerTurn, mouseClickVector)
			return framesDuration

		dt = clock.tick(fps)*timeConstant
		pygame.event.pump()

		framesDuration += 1

		# check whether the projectile needs to be despawned
		if framesDuration > framesLifeTime and projExist:
			projExist = False
			del proj

			# start the next player's turn
			playerTurn = (playerTurn+1)%numPlayers
			while not playerArray[playerTurn].playerExist:
				playerTurn = (playerTurn+1)%numPlayers

		if projExist:
			iframeCount += 1
			proj.path()

			for i in range(len(playerArray)):
				player = playerArray[i]

				# check for player collisions, and manage them accordingly
				if player.playerExist:
					# a parallelogram that represents the part
					checkPath = [prevPos[1], proj.vertices[0], proj.vertices[2], prevPos[2]]
					if collidePlayer(player.vertices, proj.vertices) or collidePlayer(player.vertices, checkPath):
						if i != playerTurn:
							player.playerExist = False
							deadPlayers += 1

						elif iframeCount > iframes:
							player.playerExist = False
							deadPlayers += 1

			# check for planet collisions and mange them accordingly
			collisionArray = collidePlanet(proj)
			if collisionArray[0]:
				phantoms.append(PhantomArrow(proj.pos, proj.v, proj.colour, collisionArray[1]))

				projExist = False
				del proj

				# set the player's turn as the next player's turn
				playerTurn = (playerTurn+1)%numPlayers
				while not playerArray[playerTurn].playerExist:
					playerTurn = (playerTurn+1)%numPlayers

		drawLoop(playerTurn, mouseClickVector)

		# check if a player is walking, and if it is legal to do so
		# update their position accordingly
		if walking and not projExist:
			playerArray[playerTurn].walk(walkDir)

		# set the position of the mouse to the middle so that it never reaches the edge of the screen
		if not mouseClick:
			pygame.mouse.set_pos(canvasWidth//2, canvasHeight//2)

		# check for and manage inputs
		for event in pygame.event.get():
			if event.type == pygame.KEYDOWN:
				if event.key in [pygame.K_a, pygame.K_d]:
					# start walking on keypress
					walkDir = event.key
					walking = True

				# esc key exits the game
				elif event.key == pygame.K_ESCAPE:
					return -1


			elif event.type == pygame.KEYUP:
				keyList = pygame.key.get_pressed()
				if not keyList[pygame.K_a] or not keyList[pygame.K_d]:
					# stop walking on keyrelease
					walking = False


			elif event.type == pygame.MOUSEBUTTONDOWN and not projExist:
				if pygame.mouse.get_pressed()[0]:
					# record the mouse position on a press of mb1
					mouseClickVector = pygame.mouse.get_pos()
					mouseClick = True

					drawForceVector = True

			elif event.type == pygame.MOUSEBUTTONUP and not projExist:
				if not pygame.mouse.get_pressed()[0] and mouseClick:
					# record the mouse position on the release of mb1
					mouseReleaseVector = pygame.mouse.get_pos()
					mouseClick = False

					# shoot a projectile accordingly
					playerArray[playerTurn].shoot(mouseClickVector, mouseReleaseVector)

					projTrail = None
					projTrail = Trail(proj.colour)

					drawForceVector = False

					# reinitialize all projectile-time variables
					iframeCount = 0
					framesDuration = 0

					prevPos = proj.pos

		if projExist:
			prevPos = proj.vertices


def endLoop(framesDuration):
	global dt
	global proj, projTrail
	global framesLifeTime

	projExist = True
	framesLifeTime //= 10

	while projExist:
		dt = clock.tick(fps)*timeConstant
		drawLoop(0,0)
		
		framesDuration += 1
		if framesDuration > framesLifeTime and projExist:
			projExist = False
			del proj

			break

		if projExist:
			proj.path()
			# check for planet collisions and mange them accordingly
			collisionArray = collidePlanet(proj)
			if collisionArray[0]:
				phantoms.append(PhantomArrow(proj.pos, proj.v, proj.colour, collisionArray[1]))

				projExist = False
				del proj

				break


if __name__ == "__main__":
	import pygame
	import os
	from random import randrange, seed
	from math import floor, pi, sin, cos, atan2
	pygame.init()

	initializeGame()
	initializePlanets()
	initializeStars()

	pygame.display.update()

	framesDuration = mainLoop()
	if framesDuration == -1:
		pygame.quit()
		quit()
	else:
		endLoop(framesDuration)

	timeAfterWin = 200

	pygame.time.delay(timeAfterWin)
	
	pygame.quit()
	quit()
