set num 200
system globSys dim $num

set count 0

while {$count < [expr $num/4]} {
	system globSys globule $count 1 1 [expr cos($count)] [expr $count * 1.1] [expr sin($count)] 0 0 0
	set count [expr $count + 1]
}
while {$count < $num/2} {
	system globSys globule $count 1 1 [expr sin($count*3.1415) * 0.5] [expr ($count - $num*0.25) * 1.1] [expr cos($count*3.1415) * 0.5] 0 0 0
	set count [expr $count + 1]
}

while {$count < [expr $num/2 + $num/4]} {
	system globSys globule $count 1 1 [expr cos($count*3.1415*0.5) * 1.4] [expr ($count - $num*0.5) * 1.1] [expr sin($count*3.1415*0.5) * 1.4] 0 0 0
	set count [expr $count + 1]
}
while {$count < $num} {
	system globSys globule $count 1 1 [expr sin($count*3.1415*0.75) * 0.5] [expr ($count - $num*0.75) * 1.1] [expr cos($count*3.1415*0.75) * 0.5] 0 0 0
	set count [expr $count + 1]
}


system globSys r0 1.0
simulator globSim constants 5.0 3.0 2.0 0.001
simulator globSim material liquid 
simulator globSim ground 5000 300
simulator globSim gravity -9.8