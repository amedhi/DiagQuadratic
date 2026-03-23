a = int(input('Enter a  number '))

print('The number entered is = ', a)

f = 1
for i in range(1, a+1):
	f = f*i
	print('i = ', i)
	print('f = ', f)
	print('\n')
