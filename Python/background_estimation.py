from numpy  import  array , uint8 , zeros
from    scipy   import  signal ,    ndimage
from VideoCapturePlayer import VideoCapturePlayer as VCP
from misc import scipyFromOpenCV

def threshold_image(np_image ,  n=[]):
    if len(n) < 5:
        # First capture a few images - give the camera time to adjust...
        n.append(np_image[:])
        return np_image

    original = n[4]
    img =   np_image [:]

    # Take the difference between the original frame and the current frame
    differenceImage = abs(np_image.astype(int) - original.astype(int)).astype(uint8)

    # amount of "Change" required before something will show up
    thresholdValue = 30
    
    # Take the N-Dimensional difference (3 channels of binary)
    differenceImage = (differenceImage >= thresholdValue)
    
    # Convert this to one channel binary
    differenceImage = differenceImage.mean(2).astype(bool)
    
    # Remove Salt & Pepper Noise
    differenceImage =   ndimage.median_filter(differenceImage , size=5)
    
    # Create a black image of same type and shape as input
    output = zeros(img.shape).astype(img.dtype)
    
    # Take the original pixel at every point where the image is "different"
    output[differenceImage] = img[differenceImage]
    return output

if __name__ == "__main__":
    title = "SciPy Background Subtraction"
    # VCP(threshold_image ,   title=title).main()
    